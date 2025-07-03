#!/bin/bash

# Source remote file copy functions if available
if [[ -f "$(dirname "$0")/remote-file-copy.sh" ]]; then
    source "$(dirname "$0")/remote-file-copy.sh"
elif [[ -f "./remote-file-copy.sh" ]]; then
    source "./remote-file-copy.sh"
fi

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
: "${EBS_QTY:="3"}" #: "normally 6"


function discover_or_launch_instances() {
    echo "Discovering existing instances and launching missing ones up to target count: ${NUM_INSTANCES}..."
    local indices_to_create=() # 0-based indices of instances to create

    for i in $(seq 0 $((NUM_INSTANCES - 1))); do
        local host_num=$((i + 1))
        local fqdn="${INSTANCE_NAME}-${host_num}.${DOMAIN}"
        local tag_name="${INSTANCE_NAME}-${host_num}"
        TARGET_FQDNS[$i]="${fqdn}"

        echo "Checking for existing instance: ${tag_name} (${fqdn})..."
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

        if [[ "${existing_instance_info}" != "null" && -n "${existing_instance_info}" ]]; then
            TARGET_INSTANCE_IDS[$i]=$(echo "${existing_instance_info}" | jq -r '.InstanceId')
            TARGET_PUBLIC_IPS[$i]=$(echo "${existing_instance_info}" | jq -r '.PublicIpAddress // ""')
            TARGET_INTERNAL_IPS[$i]=$(echo "${existing_instance_info}" | jq -r '.PrivateIpAddress // ""')
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
        echo "Need to launch ${num_to_launch} new instance(s)."

        local userdata_param=""
        if [[ -f "${CLOUD_INIT_FILE}" ]]; then
            userdata_param="--user-data file://${CLOUD_INIT_FILE}"
            echo "Using cloud-init script from ${CLOUD_INIT_FILE} for new instances."
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
            echo "Setting hostname to ${short_hostname} (${fqdn} )..."
            if ! ssh -i "${EC2_KEY_FILE}" -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null \
                "${EC2_USER}@${instance_public_ip}" \
                "sudo hostnamectl set-hostname ${short_hostname}"; then
                echo "Error: Failed to set hostname for ${short_hostname} (${instance_public_ip}). Aborting."
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
                    echo "Error: Failed to update /etc/hosts for ${internal_ip_for_host} ${fqdn} ${short_hostname} (${instance_public_ip}). Aborting."
                    exit 1
                fi
                echo "Successfully updated /etc/hosts:  ${internal_ip_for_host} ${fqdn} ${short_hostname}  -> ${internal_ip_for_host}"
            fi

            echo "Package installation (podman, lvm2) for ${fqdn} is handled by cloud-init."
        else
            echo "Skipping preparation for existing node ${TARGET_FQDNS[$i]}."
        fi
    done
    echo "All new nodes prepared successfully."
}

function wait_for_instance() {
  # This function now populates TARGET_PUBLIC_IPS and TARGET_INTERNAL_IPS for all instances in TARGET_INSTANCE_IDS
  echo "Waiting for instances to be ready..."
  local max_attempts=30
  local attempt=0
  
  while [[ ${attempt} -lt ${max_attempts} ]]; do
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
      while [[ ${ssh_attempt} -lt ${max_ssh_attempts} ]]; do
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

function bootstrap_nodes() {
    if [[ $NUM_INSTANCES -eq 0 ]]; then
        echo "No instances to bootstrap."
        return
    fi
    
    echo "Bootstrapping ${NUM_INSTANCES} nodes with bootstrap-node.sh..."
    
    # Bootstrap the first node (index 0) as the MON node
    local first_node_index=0
    local first_node_public_ip="${TARGET_PUBLIC_IPS[$first_node_index]}"
    local first_node_fqdn="${TARGET_FQDNS[$first_node_index]}"
    local first_node_internal_ip="${TARGET_INTERNAL_IPS[$first_node_index]}"
    
    if [[ -z "${first_node_public_ip}" || "${first_node_public_ip}" == "null" ]]; then
        echo "Error: First node ${first_node_fqdn} has no public IP. Cannot proceed."
        return 1
    fi
    
    echo "=== Bootstrapping first node: ${first_node_fqdn} ==="
    
    # Copy bootstrap script to the first node
    if ! scp -i "${EC2_KEY_FILE}" -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null \
        ./bootstrap-node.sh "${EC2_USER}@${first_node_public_ip}:/tmp/bootstrap-node.sh"; then
        echo "Error: Failed to copy bootstrap-node.sh to ${first_node_fqdn}."
        return 1
    fi
    
    # Execute bootstrap script on the first node with 'first' argument
    echo "Executing bootstrap-node.sh first on ${first_node_fqdn}..."
    if ! ssh -i "${EC2_KEY_FILE}" -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null \
        "${EC2_USER}@${first_node_public_ip}" "sudo HDDS_PER_SSD=${HDDS_PER_SSD:-6} bash /tmp/bootstrap-node.sh first"; then
        echo "Error: Failed to bootstrap first node ${first_node_fqdn}."
        return 1
    fi
    
    echo "First node ${first_node_fqdn} bootstrapped successfully."
    
    # Wait for SSH key to be generated
    echo "Waiting for SSH key to be available on ${first_node_fqdn}..."
    local ssh_key_attempts=30
    local ssh_key_count=0
    while ! ssh -i "${EC2_KEY_FILE}" -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null \
        "${EC2_USER}@${first_node_public_ip}" "sudo test -f /etc/ceph/ceph.pub" 2>/dev/null; do
        ssh_key_count=$((ssh_key_count + 1))
        if [[ ${ssh_key_count} -ge ${ssh_key_attempts} ]]; then
            echo "Error: SSH key did not become available after ${ssh_key_attempts} attempts."
            return 1
        fi
        echo "SSH key not yet available, waiting 10s... (Attempt ${ssh_key_count}/${ssh_key_attempts})"
        sleep 10
    done
    echo "SSH key is available on ${first_node_fqdn}."
    
    
    echo "=== Preparing additional nodes for cluster joining ==="
    
    # Prepare other nodes (they will be added to cluster via orchestrator)
    for i in $(seq 1 $((NUM_INSTANCES - 1))); do
        local target_public_ip="${TARGET_PUBLIC_IPS[$i]}"
        local target_fqdn="${TARGET_FQDNS[$i]}"
        local target_internal_ip="${TARGET_INTERNAL_IPS[$i]}"
        
        if [[ -z "${target_public_ip}" || "${target_public_ip}" == "null" ]]; then
            echo "Warning: Skipping ${target_fqdn} as its public IP is not available."
            continue
        fi
        
        echo "Running bootstrap-node.sh on ${target_fqdn} (existing or new instance)..."
        
        # Copy bootstrap script to the node
        if ! scp -i "${EC2_KEY_FILE}" -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null \
            ./bootstrap-node.sh "${EC2_USER}@${target_public_ip}:/tmp/bootstrap-node.sh"; then
            echo "Error: Failed to copy bootstrap-node.sh to ${target_fqdn}. Skipping this node."
            continue
        fi
        
        # Execute bootstrap script with 'others' argument (always run, even for existing instances)
        echo "Executing bootstrap-node.sh others on ${target_fqdn}..."
        if ! ssh -i "${EC2_KEY_FILE}" -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null \
            "${EC2_USER}@${target_public_ip}" "sudo HDDS_PER_SSD=${HDDS_PER_SSD:-6} bash /tmp/bootstrap-node.sh others"; then
            echo "Error: Failed to execute bootstrap-node.sh on ${target_fqdn}."
            continue
        fi
        
        echo "Bootstrap-node.sh completed on ${target_fqdn}."
    done
    
    echo "=== Adding nodes to cluster via orchestrator ==="
    
    # Add each additional node to cluster using the join_cluster mode
    for i in $(seq 1 $((NUM_INSTANCES - 1))); do
        local target_public_ip="${TARGET_PUBLIC_IPS[$i]}"
        local target_fqdn="${TARGET_FQDNS[$i]}"
        local target_internal_ip="${TARGET_INTERNAL_IPS[$i]}"
        local target_hostname="${target_fqdn%%.*}"  # Extract short hostname
        
        if [[ -z "${target_public_ip}" || "${target_public_ip}" == "null" ]]; then
            echo "Warning: Skipping ${target_fqdn} as its public IP is not available."
            continue
        fi
        
        echo "Adding ${target_fqdn} to cluster via join_cluster mode..."
        
        # Copy bootstrap script to the node
        if ! scp -i "${EC2_KEY_FILE}" -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null \
            ./bootstrap-node.sh "${EC2_USER}@${target_public_ip}:/tmp/bootstrap-node.sh"; then
            echo "Error: Failed to copy bootstrap-node.sh to ${target_fqdn}. Skipping this node."
            continue
        fi
        
        # Use join_cluster mode to add the node to cluster from the first node
        echo "Executing join_cluster mode on first node for ${target_fqdn}..."
        if ssh -i "${EC2_KEY_FILE}" -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null \
            "${EC2_USER}@${first_node_public_ip}" \
            "sudo FIRST_NODE_INTERNAL_IP=${first_node_internal_ip} TARGET_HOSTNAME=${target_hostname} TARGET_INTERNAL_IP=${target_internal_ip} EC2_USER=${EC2_USER} FIRST_NODE_PUBLIC_IP=${first_node_public_ip} TARGET_PUBLIC_IP=${target_public_ip} EC2_KEY_FILE=${EC2_KEY_FILE} bash /tmp/bootstrap-node.sh join_cluster"; then
            echo "Successfully added ${target_fqdn} to cluster."
            
            # Wait for cluster configuration to be distributed to the new node  
            echo "Waiting for cluster configuration to be available on ${target_fqdn}..."
            local config_attempts=30
            local config_count=0
            while [[ ${config_count} -lt ${config_attempts} ]]; do
                if ssh -i "${EC2_KEY_FILE}" -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null \
                    "${EC2_USER}@${target_public_ip}" "sudo test -f /etc/ceph/ceph.conf" 2>/dev/null; then
                    echo "Cluster configuration is now available on ${target_fqdn}."
                    break
                fi
                
                config_count=$((config_count + 1))
                echo "Waiting for cluster config on ${target_fqdn}... (Attempt ${config_count}/${config_attempts})"
                sleep 10
            done
            
            if [[ ${config_count} -eq ${config_attempts} ]]; then
                echo "Warning: Cluster configuration did not become available on ${target_fqdn} within expected time."
                echo "Skipping OSD creation for this node."
            else
                # Now create OSDs on this node
                echo "Creating OSDs on ${target_fqdn} (cluster config available)..."
                
                # Execute bootstrap script with 'others' argument to create OSDs
                if ssh -i "${EC2_KEY_FILE}" -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null \
                    "${EC2_USER}@${target_public_ip}" "sudo HDDS_PER_SSD=${HDDS_PER_SSD:-6} bash /tmp/bootstrap-node.sh others"; then
                    echo "OSDs created successfully on ${target_fqdn}."
                else
                    echo "Warning: Failed to create OSDs on ${target_fqdn}."
                fi
            fi
        else
            echo "Warning: Failed to add ${target_fqdn} to cluster via join_cluster mode."
            echo "Manual addition may be required: ceph orch host add ${target_hostname} ${target_internal_ip}"
        fi
    done
    
    # Clean up bootstrap scripts on all nodes
    for i in $(seq 0 $((NUM_INSTANCES - 1))); do
        local target_public_ip="${TARGET_PUBLIC_IPS[$i]}"
        if [[ -n "${target_public_ip}" && "${target_public_ip}" != "null" ]]; then
            ssh -i "${EC2_KEY_FILE}" -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null \
                "${EC2_USER}@${target_public_ip}" "rm -f /tmp/bootstrap-node.sh" 2>/dev/null || true
        fi
    done
    
    echo "=== OSD creation completed inline after each node joined cluster ==="
    
    echo "Multi-node cluster deployment complete."
}

# Attempt to retrieve AWS identity information
identity_info=$(aws sts get-caller-identity --query '[Account, Arn]' --output text 2>/dev/null)
# Check if the command succeeded
if [[ $? -ne 0 ]]; then
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
bootstrap_nodes

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

