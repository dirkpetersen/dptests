#!/bin/bash

# AWS EC2 Instance Launch and Configuration Script
#
# This script automates the process of launching one or more AWS EC2 instances.
# Key functionalities include:
# - Launching a specified number of EC2 instances with configurable parameters
#   (region, instance type, AMI, root volume size, instance name, domain, etc.).
# - Utilizing a cloud-init script (if provided) for initial instance setup.
# - Checking AWS CLI authentication and attempting SSO login if needed.
# - Managing EC2 key pairs for SSH access.
# - Configuring storage:
#   - Setting the root volume size.
#   - Attaching a specified number of additional EBS volumes to each instance
#     by calling the 'ebs-create-attach.sh' script.
# - Waiting for instances to reach a 'running' state and for SSH to become available.
# - Optionally registering DNS A records for each instance in AWS Route 53
#   if a hosted zone is found.
# - Outputting SSH commands to connect to the newly created instances.
#
# Usage:
#   ./ec2-create-instances.sh [number_of_instances]
#
# Environment variables can be used to override default configurations, e.g.:
#   AWS_REGION, EC2_TYPE, AMI_IMAGE, ROOT_VOLUME_SIZE, INSTANCE_NAME,
#   DOMAIN, CLOUD_INIT_FILE, EC2_USER, EC2_SECURITY_GROUP, EBS_TYPE,
#   EBS_SIZE, EBS_QTY, EC2_KEY_NAME, EC2_KEY_FILE, AWS_AZ.

declare -a TARGET_INSTANCE_IDS # Array to store instance IDs for all target slots
declare -a TARGET_PUBLIC_IPS   # Array to store public IPs for all target slots
declare -a TARGET_FQDNS        # Array to store FQDNs for all target slots
declare -a IS_INSTANCE_NEW     # Boolean array: true if instance at index was created in this run
declare -a TARGET_INTERNAL_IPS # Array to store internal IPs for all target slots

# AWS EC2 Instance Launch Script
# This script launches a c7gd.medium EC2 instance with Rocky Linux 9

# Get number of instances from command line
NUM_INSTANCES=${1:-1}

if [[ ! "$NUM_INSTANCES" =~ ^[1-9][0-9]*$ ]]; then
    echo "Error: Number of instances must be a positive integer"
    exit 1
fi

# Configuration Variables
: "${AWS_REGION:="us-west-2"}"
: "${EC2_TYPE:="c7gd.medium"}"
: "${AMI_IMAGE:="ami-03be04a3da3a40226"}"
: "${ROOT_VOLUME_SIZE:="16"}"
: "${INSTANCE_NAME:="ceph-test"}"
: "${DOMAIN:="ai.oregonstate.edu"}"
: "${CLOUD_INIT_FILE:="ec2-cloud-init.txt"}"
: "${EC2_USER:="rocky"}"
: "${EC2_SECURITY_GROUP:="SSH-HTTP-ICMP"}"
: "${EBS_TYPE:="st1"}"
: "${EBS_SIZE:="125"}"
: "${EBS_QTY:="2"}" #: "normally 6"


function discover_or_launch_instances() {
    echo "Discovering existing instances and launching missing ones up to target count: $NUM_INSTANCES..."
    local indices_to_create=() # 0-based indices of instances to create

    for i in $(seq 0 $((NUM_INSTANCES - 1))); do
        local host_num=$((i + 1))
        local fqdn="${INSTANCE_NAME}-${host_num}.${DOMAIN}"
        local tag_name="${INSTANCE_NAME}-${host_num}"
        TARGET_FQDNS[$i]="$fqdn"

        echo "Checking for existing instance: $tag_name ($fqdn)..."
        local existing_instance_info
        local aws_cli_output
        local aws_cli_stderr_file
        aws_cli_stderr_file=$(mktemp)

        if ! aws_cli_output=$(aws ec2 describe-instances \
            --filters "Name=tag:Name,Values=${tag_name}" "Name=instance-state-name,Values=pending,running" \
            --query "Reservations[].Instances[] | [?InstanceLifecycle!='spot'].{InstanceId:InstanceId,PublicIpAddress:PublicIpAddress,PrivateIpAddress:PrivateIpAddress,LaunchTime:LaunchTime}" \
            --output json 2> "${aws_cli_stderr_file}"); then
            echo "Warning: AWS CLI command failed for tag ${tag_name}. Error: $(cat "${aws_cli_stderr_file}")"
            existing_instance_info="null" # Treat as not found
        else
            # If AWS CLI succeeded, aws_cli_output contains the JSON
            local jq_stderr_file
            jq_stderr_file=$(mktemp)
            if ! existing_instance_info=$(echo "${aws_cli_output}" | jq -c 'sort_by(.LaunchTime) | .[0] // null' 2> "${jq_stderr_file}"); then
                echo "Warning: jq command failed for tag ${tag_name}. Error: $(cat "${jq_stderr_file}")"
                echo "DEBUG: Input to jq for tag ${tag_name} was: ${aws_cli_output}" # Print the input that made jq fail
                existing_instance_info="null" # Treat as not found
            fi
            rm -f "${jq_stderr_file}"
        fi
        rm -f "${aws_cli_stderr_file}"

        if [[ "$existing_instance_info" != "null" && -n "$existing_instance_info" ]]; then
            TARGET_INSTANCE_IDS[$i]=$(echo "$existing_instance_info" | jq -r '.InstanceId')
            TARGET_PUBLIC_IPS[$i]=$(echo "$existing_instance_info" | jq -r '.PublicIpAddress // ""')
            TARGET_INTERNAL_IPS[$i]=$(echo "$existing_instance_info" | jq -r '.PrivateIpAddress // ""')
            IS_INSTANCE_NEW[$i]=false
            echo "Found existing instance ${tag_name} (ID: ${TARGET_INSTANCE_IDS[$i]})."
        else
            echo "Instance ${tag_name} not found or not in a usable state. Will be created."
            indices_to_create+=($i)
            IS_INSTANCE_NEW[$i]=true
            TARGET_INSTANCE_IDS[$i]="" # Placeholder
            TARGET_PUBLIC_IPS[$i]="" # Placeholder
            TARGET_INTERNAL_IPS[$i]="" # Placeholder
        fi
    done

    # Determine and configure the security group REGARDLESS of whether new instances are launched.
    # This security group will be used by all instances, new or existing.
    echo "Fetching security group ID for: ${EC2_SECURITY_GROUP}"
    local security_group_id # Declare here so it's available for launch if needed
    security_group_id=$(aws ec2 describe-security-groups \
        --region "${AWS_REGION}" \
        --filters "Name=group-name,Values=${EC2_SECURITY_GROUP}" \
        --query 'SecurityGroups[0].GroupId' \
        --output text)

    if [[ -z "$security_group_id" ]]; then
        echo "Error: Security group '${EC2_SECURITY_GROUP}' not found in region '${AWS_REGION}'."
        exit 1
    fi
    echo "Using security group ID: ${security_group_id} (Name: ${EC2_SECURITY_GROUP})"

    # Ensure the security group allows SSH from itself for intra-cluster communication (e.g., Ceph orchestration)
    # Check if the rule already exists: TCP, port 22, source is the security group itself.
    local rule_exists
    rule_exists=$(aws ec2 describe-security-groups \
        --region "${AWS_REGION}" \
        --group-ids "${security_group_id}" \
        --filters "Name=ip-permission.protocol,Values=tcp" \
                  "Name=ip-permission.from-port,Values=22" \
                  "Name=ip-permission.to-port,Values=22" \
                  "Name=ip-permission.user-id-group-pairs.group-id,Values=${security_group_id}" \
        --query "SecurityGroups[0].IpPermissions[?FromPort==\`22\` && ToPort==\`22\` && IpProtocol=='tcp' && UserIdGroupPairs[?GroupId=='${security_group_id}']].UserIdGroupPairs" \
        --output text 2>/dev/null)

    if [[ -z "$rule_exists" ]]; then
        echo "Adding ingress rule to security group ${security_group_id} to allow SSH (port 22) from itself..."
        local auth_stderr_file
        auth_stderr_file=$(mktemp)
        if aws ec2 authorize-security-group-ingress \
            --region "${AWS_REGION}" \
            --group-id "${security_group_id}" \
            --protocol tcp \
            --port 22 \
            --source-group "${security_group_id}" 2> "${auth_stderr_file}"; then
            echo "Successfully added SSH ingress rule (port 22, source: self) to security group ${security_group_id}."
        else
            # Check if the error was because the rule already exists
            if grep -q "InvalidPermission.Duplicate" "${auth_stderr_file}"; then
                echo "SSH ingress rule (port 22, source: self) already exists in security group ${security_group_id} (confirmed by authorize attempt; this is not an error)."
            else
                # A real error occurred
                echo "Error: Failed to add SSH ingress rule to security group ${security_group_id}. Ceph internode SSH might fail."
                echo "AWS CLI Error: $(cat "${auth_stderr_file}")"
                echo "Please manually add an Inbound rule: TCP, Port 22, Source: ${security_group_id}"
                # Consider exiting here if this rule is critical and automation fails:
                # exit 1
            fi
        fi
        rm -f "${auth_stderr_file}"
    else
        echo "SSH ingress rule (port 22, source: self) already exists in security group ${security_group_id}."
    fi

    local num_to_launch=${#indices_to_create[@]}
    if [[ $num_to_launch -gt 0 ]]; then
        echo "Need to launch $num_to_launch new instance(s)."

        local userdata_param=""
        if [[ -f "$CLOUD_INIT_FILE" ]]; then
            userdata_param="--user-data file://${CLOUD_INIT_FILE}"
            echo "Using cloud-init script from $CLOUD_INIT_FILE for new instances."
        fi

        local blk_dev_json="["
        if [[ $ROOT_VOLUME_SIZE -gt 0 ]]; then
            local root_dev_name=$(aws ec2 describe-images --image-ids ${AMI_IMAGE} --query 'Images[].RootDeviceName' --output text)
            blk_dev_json+="{\"DeviceName\":\"${root_dev_name}\",\"Ebs\":{\"VolumeSize\":${ROOT_VOLUME_SIZE}}}"
            # Assuming ephemeral0 is desired if root volume is customized. Adjust if not.
            blk_dev_json+=",{\"DeviceName\":\"/dev/sdb\",\"VirtualName\":\"ephemeral0\"}"
        fi
        blk_dev_json+="]"

        local az_param=""
        if [[ -n "${AWS_AZ}" ]]; then
            az_param="--placement AvailabilityZone=${AWS_AZ}"
            echo "Using specified availability zone for new instances: ${AWS_AZ}"
        fi

        echo "Launching $num_to_launch instances..."
        local launch_output
        launch_output=$(aws ec2 run-instances --region ${AWS_REGION} \
            --image-id ${AMI_IMAGE} \
            --count ${num_to_launch} \
            --instance-type ${EC2_TYPE} \
            --key-name ${EC2_KEY_NAME} \
            --security-group-ids "${security_group_id}" \
            --block-device-mappings "${blk_dev_json}" \
            ${userdata_param} \
            ${az_param} \
            --output json)

        if [[ -z "$launch_output" ]]; then
            echo "Error: Failed to launch new instances. Check AWS credentials and permissions."
            exit 1
        fi

        local launched_instance_ids_from_aws=($(echo "$launch_output" | jq -r '.Instances[].InstanceId'))
        if [[ ${#launched_instance_ids_from_aws[@]} -ne $num_to_launch ]]; then
            echo "Error: Expected $num_to_launch new instance IDs, but got ${#launched_instance_ids_from_aws[@]}."
            exit 1
        fi

        for j in $(seq 0 $((num_to_launch - 1))); do
            local target_array_index=${indices_to_create[$j]}
            local new_instance_id=${launched_instance_ids_from_aws[$j]}
            TARGET_INSTANCE_IDS[$target_array_index]=$new_instance_id
            local tag_name="${INSTANCE_NAME}-$((target_array_index + 1))"
            aws ec2 create-tags --resources ${new_instance_id} --tags "Key=Name,Value=${tag_name}"
            echo "Launched new instance ${tag_name} (ID: ${new_instance_id}). It will be fully configured."
        done
    else
        echo "All $NUM_INSTANCES target instances already exist."
    fi
    echo "Instance discovery and launch phase complete."
}

function add_disks() {
    echo "Adding disks to newly created instances..."
    for i in $(seq 0 $((NUM_INSTANCES - 1))); do
        if [[ "${IS_INSTANCE_NEW[$i]}" == true ]]; then
            local instance_id=${TARGET_INSTANCE_IDS[$i]}
            echo "Adding ${EBS_QTY} ${EBS_TYPE} volumes (${EBS_SIZE}GB) to new instance ${TARGET_FQDNS[$i]} (ID: ${instance_id})"
            ./ebs-create-attach.sh "${instance_id}" "${EBS_TYPE}" "${EBS_SIZE}" "${EBS_QTY}"
        fi
    done
}

function register_dns() {
  hosted_zone_id=$(aws route53 list-hosted-zones --query 'HostedZones[0].Id' --output text)
  if [[ -n "$hosted_zone_id" && "$hosted_zone_id" != "None" ]]; then
    for i in $(seq 0 $((NUM_INSTANCES - 1))); do
      local fqdn="${TARGET_FQDNS[$i]}"
      local public_ip="${TARGET_PUBLIC_IPS[$i]}"

      if [[ -z "$public_ip" || "$public_ip" == "null" ]]; then
          echo "Skipping DNS registration for ${fqdn} as public IP is not available."
          continue
      fi
      
      change_batch=$(cat <<EOF
{
  "Changes": [
    {
      "Action": "UPSERT",
      "ResourceRecordSet": {
        "Name": "${fqdn}",
        "Type": "A",
        "TTL": 300,
        "ResourceRecords": [
          {
            "Value": "${public_ip}"
          }
        ]
      }
    }
  ]
}
EOF
      )

      if aws route53 change-resource-record-sets \
        --hosted-zone-id $hosted_zone_id \
        --change-batch "$change_batch" > /dev/null; then
        echo "DNS record UPSERTED for ${fqdn} -> ${public_ip}"
      else
        echo "Failed to UPSERT DNS record for ${fqdn} -> ${public_ip}"
      fi
    done
  else
    echo "No Route 53 hosted zone found. Using public IP addresses."
  fi
}

function prepare_new_nodes() {
    echo "Preparing all new instances (setting hostname, installing packages)..."
    for i in $(seq 0 $((NUM_INSTANCES - 1))); do
        if [[ "${IS_INSTANCE_NEW[$i]}" == true ]]; then
            local instance_public_ip="${TARGET_PUBLIC_IPS[$i]}"
            local fqdn="${TARGET_FQDNS[$i]}"
            
            echo "Preparing new node ${fqdn} (${instance_public_ip})..."

            # Extract short hostname from FQDN
            local short_hostname="${fqdn%%.*}"
            
            # Set hostname
            echo "Setting hostname to ${fqdn}..."
            if ! ssh -i "${EC2_KEY_FILE}" -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null \
                "${EC2_USER}@${instance_public_ip}" \
                "sudo hostnamectl set-hostname ${fqdn}"; then
                echo "Error: Failed to set hostname for ${fqdn} (${instance_public_ip}). Aborting."
                exit 1 # Exit script if hostname setting fails
            fi
            echo "Successfully set hostname for ${fqdn} (short name: ${short_hostname})."

            # Update /etc/hosts to ensure FQDN resolves to external IP and short hostname to internal IP
            local internal_ip_for_host="${TARGET_INTERNAL_IPS[$i]}"
            if [[ -n "$internal_ip_for_host" && "$internal_ip_for_host" != "null" ]]; then
                echo "Updating /etc/hosts for ${fqdn} and ${short_hostname}..."
                # do not use: echo '${instance_public_ip} ${fqdn}' | sudo tee -a /etc/hosts; \
                # Remove old entries for this FQDN or short hostname, then add the new ones                
                local update_hosts_cmd="sudo sed -i -e '/\\s${fqdn}\$/d' -e '/\\s${short_hostname}\$/d' /etc/hosts; \                
                echo '${internal_ip_for_host} ${fqdn} ${short_hostname}' | sudo tee -a /etc/hosts"
                if ! ssh -i "${EC2_KEY_FILE}" -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null \
                    "${EC2_USER}@${instance_public_ip}" "${update_hosts_cmd}"; then
                    echo "Error: Failed to update /etc/hosts for ${fqdn} (${instance_public_ip}). Aborting."
                    exit 1
                fi
                echo "Successfully updated /etc/hosts: ${fqdn} -> ${instance_public_ip}, ${short_hostname} -> ${internal_ip_for_host}"
            fi

            echo "Package installation (podman, lvm2) for ${fqdn} is handled by cloud-init."
        else
            echo "Skipping preparation for existing node ${TARGET_FQDNS[$i]}."
        fi
    done
    echo "All new nodes prepared successfully."
}

function configure_ceph_nodes() {
    if [[ $NUM_INSTANCES -eq 0 ]]; then
        echo "No target instances to configure for Ceph."
        return
    fi

    local key_source_node_ip=""
    local key_source_node_fqdn_display="" # For logging
    local temp_ceph_pub_key_local="/tmp/ceph_cluster_key_$(date +%s%N).pub" # Unique temp file on control machine
    
    # The first instance (index 0) is always the MON/bootstrap/key source node.
    local mon_node_index=0
    key_source_node_ip="${TARGET_PUBLIC_IPS[$mon_node_index]}"
    key_source_node_fqdn_display="${TARGET_FQDNS[$mon_node_index]}"

    if [[ -z "$key_source_node_ip" || "$key_source_node_ip" == "null" ]]; then
        echo "Error: MON node ${key_source_node_fqdn_display} does not have a public IP. Cannot proceed with Ceph configuration."
        return 1
    fi

    if [[ "${IS_INSTANCE_NEW[$mon_node_index]}" == true ]]; then
        echo "MON node ${key_source_node_fqdn_display} is newly created. Bootstrapping Ceph..."

        echo "Copying ./ceph-bootstrap.sh to ${EC2_USER}@${key_source_node_ip}:/tmp/ceph-bootstrap.sh..."
        scp -i "${EC2_KEY_FILE}" -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null \
            ./ceph-bootstrap.sh "${EC2_USER}@${key_source_node_ip}:/tmp/ceph-bootstrap.sh"
        if [[ $? -ne 0 ]]; then
            echo "Error: Failed to copy ceph-bootstrap.sh to $key_source_node_fqdn_display. Aborting Ceph configuration."
            return 1
        fi

        # Wait for podman to be installed by cloud-init
        echo "Waiting for podman to be available on ${key_source_node_fqdn_display} (installed by cloud-init)..."
        local podman_check_attempts=30 # Approx 5 minutes (30 * 10s)
        local podman_check_count=0
        while ! ssh -i "${EC2_KEY_FILE}" -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null \
                "${EC2_USER}@${key_source_node_ip}" "command -v podman" >/dev/null 2>&1; do
            podman_check_count=$((podman_check_count + 1))
            if [[ $podman_check_count -ge $podman_check_attempts ]]; then
                echo "Error: podman did not become available on ${key_source_node_fqdn_display} after ${podman_check_attempts} attempts."
                return 1
            fi
            echo "podman not yet available, waiting 10s... (Attempt ${podman_check_count}/${podman_check_attempts})"
            sleep 10
        done
        echo "podman is available on ${key_source_node_fqdn_display}."

        echo "Running ceph-bootstrap.sh on $key_source_node_fqdn_display..."
        ssh -A -i "${EC2_KEY_FILE}" -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null \
            "${EC2_USER}@${key_source_node_ip}" \
            "sudo bash /tmp/ceph-bootstrap.sh"
        if [[ $? -ne 0 ]]; then
            echo "Error: Failed to execute ceph-bootstrap.sh on $key_source_node_fqdn_display. Aborting Ceph configuration."
            # Attempt to clean up bootstrap script on failure
            ssh -i "${EC2_KEY_FILE}" -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null \
                "${EC2_USER}@${key_source_node_ip}" "rm -f /tmp/ceph-bootstrap.sh"
            return 1
        fi
    else
        echo "MON node ${key_source_node_fqdn_display} already exists. Using it as key source for Ceph public key."
    fi

    # === Fetch Ceph Public Key from MON node ===
    echo "Fetching Ceph public key from $key_source_node_fqdn_display ($key_source_node_ip)..."

    # Prepare Ceph public key on the key_source_node (copy to /tmp, chown, chmod)
    # This requires SSH access to key_source_node_ip using EC2_KEY_FILE (or appropriate auth for existing admin)
    echo "Preparing Ceph public key on $key_source_node_fqdn_display..."
    local prep_key_cmd="sudo cp /etc/ceph/ceph.pub /tmp/ceph.pub && sudo chown ${EC2_USER}:${EC2_USER} /tmp/ceph.pub && sudo chmod 644 /tmp/ceph.pub"
    ssh -i "${EC2_KEY_FILE}" -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null \
        "${EC2_USER}@${key_source_node_ip}" "${prep_key_cmd}"
    if [[ $? -ne 0 ]]; then
        echo "Error: Failed to prepare Ceph public key on $key_source_node_fqdn_display (/tmp/ceph.pub). Aborting key distribution."
        if [[ "${IS_INSTANCE_NEW[$mon_node_index]}" == true ]]; then # If MON was new, cleanup bootstrap script
            ssh -i "${EC2_KEY_FILE}" -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null \
                "${EC2_USER}@${TARGET_PUBLIC_IPS[$mon_node_index]}" "rm -f /tmp/ceph-bootstrap.sh"
        fi
        return 1
    fi

    # Copy Ceph public key from the key_source_node to the control machine
    echo "Fetching Ceph public key from $key_source_node_fqdn_display to control machine..."
    scp -i "${EC2_KEY_FILE}" -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null \
        "${EC2_USER}@${key_source_node_ip}:/tmp/ceph.pub" "${temp_ceph_pub_key_local}"
    if [[ $? -ne 0 ]]; then
        echo "Error: Failed to copy Ceph public key from $key_source_node_fqdn_display to control machine. Aborting key distribution."
        # Cleanup on key_source_node and potentially bootstrap script
        ssh -i "${EC2_KEY_FILE}" -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null \
            "${EC2_USER}@${key_source_node_ip}" "rm -f /tmp/ceph.pub"
        if [[ "${IS_INSTANCE_NEW[$mon_node_index]}" == true ]]; then
            ssh -i "${EC2_KEY_FILE}" -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null \
                "${EC2_USER}@${TARGET_PUBLIC_IPS[$mon_node_index]}" "rm -f /tmp/ceph-bootstrap.sh"
        fi
        return 1
    fi

    # Clean up /tmp/ceph.pub on the key_source_node as it's now on the control machine
    ssh -i "${EC2_KEY_FILE}" -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null \
        "${EC2_USER}@${key_source_node_ip}" "rm -f /tmp/ceph.pub"

    # Distribute Ceph public key to ALL target nodes (including MON itself, ssh-copy-id is idempotent)
    # This ensures all nodes (new or existing that are part of the target count) have the MON's Ceph key.
    if [[ $NUM_INSTANCES -gt 0 ]]; then
        echo "Distributing Ceph public key from ${key_source_node_fqdn_display} to all ${NUM_INSTANCES} target nodes..."
        for i in $(seq 0 $((NUM_INSTANCES - 1))); do
            local target_public_ip="${TARGET_PUBLIC_IPS[$i]}"
            local target_fqdn_display="${TARGET_FQDNS[$i]}"

            if [[ -z "$target_public_ip" || "$target_public_ip" == "null" ]]; then
                echo "Warning: Skipping key distribution to ${target_fqdn_display} as its public IP is not available."
                continue
            fi
            # Skip distributing key to the MON node itself if it was the source and already exists (no need to ssh-copy-id to self)
            # However, ssh-copy-id to self is harmless and ensures the key is there if somehow removed.
            # For simplicity and robustness, we'll attempt to copy to all, including the MON.
            echo "Distributing Ceph public key to ${target_fqdn_display} (${target_public_ip})..."
            
            # ssh-copy-id from control machine directly to target_public_ip
            # EC2_KEY_FILE is used for authentication.
            local ssh_options_for_copy_id=(
                "-o" "StrictHostKeyChecking=no"
                "-o" "UserKnownHostsFile=/dev/null"
                "-o" "IdentityFile=${EC2_KEY_FILE}"
                "-o" "ConnectTimeout=10"
            )

            if ssh-copy-id -f -i "${temp_ceph_pub_key_local}" "${ssh_options_for_copy_id[@]}" "${EC2_USER}@${target_public_ip}"; then
                echo "Successfully copied Ceph public key to ${target_fqdn_display} (${target_public_ip})."
            else
                echo "Warning: Failed to copy Ceph public key to ${target_fqdn_display} (${target_public_ip}). Exit status: $?"
            fi
        done
    else
        echo "No nodes to distribute Ceph public key to."
    fi
    
    # Final Cleanup
    rm -f "${temp_ceph_pub_key_local}" # Remove temp key from control machine
    if [[ "${IS_INSTANCE_NEW[$mon_node_index]}" == true ]]; then # If MON was new, cleanup bootstrap script
        ssh -i "${EC2_KEY_FILE}" -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null \
            "${EC2_USER}@${TARGET_PUBLIC_IPS[$mon_node_index]}" "rm -f /tmp/ceph-bootstrap.sh"
    fi
    echo "Ceph node configuration process complete."
}

function wait_for_instance() {
  # This function now populates TARGET_PUBLIC_IPS and TARGET_INTERNAL_IPS for all instances in TARGET_INSTANCE_IDS
  echo "Waiting for instances to be ready..."
  local max_attempts=30
  local attempt=0
  
  while [ $attempt -lt $max_attempts ]; do
      all_ready=true
      for i in $(seq 0 $((NUM_INSTANCES - 1))); do
          local instance_id=${TARGET_INSTANCE_IDS[$i]}
          if [[ -z "$instance_id" ]]; then # Should not happen if discover_or_launch_instances worked
              all_ready=false; break
          fi
          
          local instance_info
          instance_info=$(aws ec2 describe-instances --region ${AWS_REGION} --instance-ids "$instance_id" \
              --query 'Reservations[0].Instances[0].{State:State.Name,PublicIpAddress:PublicIpAddress,PrivateIpAddress:PrivateIpAddress,InstanceType:InstanceType}' \
              --output json 2>/dev/null)
          
          if [[ $? -eq 0 ]]; then
              local state=$(echo "$instance_info" | jq -r '.State')
              local public_ip=$(echo "$instance_info" | jq -r '.PublicIpAddress // ""')
              local private_ip=$(echo "$instance_info" | jq -r '.PrivateIpAddress // ""')
              TARGET_PUBLIC_IPS[$i]="$public_ip" # Update public IP
              TARGET_INTERNAL_IPS[$i]="$private_ip" # Update private IP
              
              if [[ "$state" != "running" || "$public_ip" == "null" || -z "$public_ip" ]]; then
                  all_ready=false
                  break
              fi
          fi
      done
      
      if $all_ready; then
          echo "All instances are running"
          break
      fi
      
      echo "Waiting for instances to be ready... (Attempt $((attempt+1))/$max_attempts)"
      sleep 10
      attempt=$((attempt+1))
  done

  if [[ $attempt -eq $max_attempts ]]; then
      echo "Error: Not all instances became ready within the expected time."
      exit 1
  fi

  echo "Waiting for SSH to become available on all instances..."
  for i in $(seq 0 $((NUM_INSTANCES - 1))); do
      local public_ip=${TARGET_PUBLIC_IPS[$i]}
      local fqdn=${TARGET_FQDNS[$i]}
      local instance_id=${TARGET_INSTANCE_IDS[$i]}

      if [[ -z "$public_ip" || "$public_ip" == "null" ]]; then
          # If an existing instance lost its public IP or a new one failed to get one.
          echo "Warning: Public IP for ${fqdn} (ID: ${instance_id}) is not available. SSH check might fail or be skipped."
          # If it's a new instance, this is a critical failure.
          if [[ "${IS_INSTANCE_NEW[$i]}" == true ]]; then
              echo "Error: Newly created instance ${fqdn} (ID: ${instance_id}) has no public IP. Aborting."
              exit 1
          fi
          continue # Skip SSH check for existing instances without public IP if that's acceptable.
      fi

      max_ssh_attempts=30
      ssh_attempt=0
      sshout=$(mktemp -t ec2-XXX)
      while [ $ssh_attempt -lt $max_ssh_attempts ]; do
          if ssh -i "${EC2_KEY_FILE}" -o ConnectTimeout=5 -o BatchMode=yes -o StrictHostKeyChecking=no ${EC2_USER}@${public_ip} exit 2>${sshout}; then
              echo "SSH is available on $public_ip"
              break
          fi
          printf "Debug SSH to ${fqdn} ($public_ip): " && cat $sshout
          echo "Waiting for SSH on ${fqdn}... (Attempt $((ssh_attempt+1))/$max_ssh_attempts)"
          sleep 10
          ssh_attempt=$((ssh_attempt+1))
      done

      if [[ $ssh_attempt -eq $max_ssh_attempts ]]; then
          echo "Error: SSH did not become available for ${fqdn} ($public_ip) within the expected time."
          exit 1
      fi
  done
}

# Attempt to retrieve AWS identity information
identity_info=$(aws sts get-caller-identity --query '[Account, Arn]' --output text 2>/dev/null)
# Check if the command succeeded
if [ $? -ne 0 ]; then
    echo "Error: AWS CLI is not authenticated. Executing 'aws sso login --no-browser' or set up your credentials."
    aws sso login --no-browser
    exit 1
fi
# Extract the account ID and the username/role name if the command succeeded
AWSACCOUNT=$(echo "$identity_info" | awk '{print $1}')
AWSUSER=$(echo "$identity_info" | awk '{print $2}' | awk -F'/' '{print $NF}')
AWSUSER2="${AWSUSER%@*}"
  # Set EC2_KEY_NAME and EC2_KEY_FILE
: "${EC2_KEY_NAME:="auto-ec2-${AWSUSER}"}"
: "${EC2_KEY_FILE:="~/.ssh/auto-ec2-${AWSACCOUNT}-${AWSUSER2}.pem"}"
EC2_KEY_FILE=$(eval echo "${EC2_KEY_FILE}")

discover_or_launch_instances
wait_for_instance
prepare_new_nodes
add_disks
register_dns
configure_ceph_nodes

echo -e "\nInstances created on ${EC2_TYPE} with AMI ${AMI_IMAGE}"

echo -e "\nInstance Information:"
for i in $(seq 0 $((NUM_INSTANCES - 1))); do
    fqdn="${TARGET_FQDNS[$i]}"
    echo "${fqdn}, external IP: ${TARGET_PUBLIC_IPS[$i]:-N/A}, internal IP: ${TARGET_INTERNAL_IPS[$i]:-N/A}"
done

echo -e "\nSSH commands to connect:"
FILE2=~${EC2_KEY_FILE#"$HOME"}
for i in $(seq 0 $((NUM_INSTANCES - 1))); do
    fqdn="${TARGET_FQDNS[$i]}"
    echo "ssh -i '${FILE2}' ${EC2_USER}@${fqdn}"
done

