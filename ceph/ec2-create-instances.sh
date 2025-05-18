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

declare -a public_ips

# AWS EC2 Instance Launch Script
# This script launches a c7gd.medium EC2 instance with Rocky Linux 9

# Get number of instances from command line
NUM_INSTANCES=${1:-1}
EXISTING_CEPH_ADMIN_FQDN=${2:-""} # Second argument for existing Ceph admin node FQDN

if [[ ! "$NUM_INSTANCES" =~ ^[1-9][0-9]*$ ]]; then
    echo "Error: Number of instances must be a positive integer"
    exit 1
fi
if [[ -n "$EXISTING_CEPH_ADMIN_FQDN" ]]; then
    echo "Cluster expansion mode: Will add Ceph public key from $EXISTING_CEPH_ADMIN_FQDN to $NUM_INSTANCES new node(s)."
fi

# Configuration Variables
: "${AWS_REGION:="us-west-2"}"
: "${EC2_TYPE:="c7gd.medium"}"
: "${AMI_IMAGE:="ami-03be04a3da3a40226"}"
: "${ROOT_VOLUME_SIZE:="16"}"
: "${INSTANCE_NAME:="ceph-test"}"
: "${DOMAIN:="ai.oregonstate.edu"}"
: "${CLOUD_INIT_FILE:="ec2-cloud-init.txt.1"}"
: "${EC2_USER:="rocky"}"
: "${EC2_SECURITY_GROUP:="SSH-HTTP-ICMP"}"
: "${EBS_TYPE:="st1"}"
: "${EBS_SIZE:="125"}"
: "${EBS_QTY:="6"}"

function launch_instance() {
  # Check for ec2-cloud-init.txt file
  userdata=""
  if [[ -f "$CLOUD_INIT_FILE" ]]; then  
    userdata="--user-data file://${CLOUD_INIT_FILE}"
    echo "Using cloud-init script from $CLOUD_INIT_FILE"
  fi

  # Construct block device mapping JSON
  BLK_DEV_JSON="["
  if [[ $ROOT_VOLUME_SIZE -gt 0 ]]; then
      ROOT_DEV=$(aws ec2 describe-images --image-ids ${AMI_IMAGE} --query 'Images[].RootDeviceName' --output text)
      BLK_DEV_JSON+="{\"DeviceName\":\"${ROOT_DEV}\",\"Ebs\":{\"VolumeSize\":${ROOT_VOLUME_SIZE}}}"
      [[ $ROOT_VOLUME_SIZE -gt 0 ]] && BLK_DEV_JSON+=","
      BLK_DEV_JSON+="{\"DeviceName\":\"/dev/sdb\",\"VirtualName\":\"ephemeral0\"}"
  fi
  BLK_DEV_JSON+="]"

  # Add availability zone parameter if specified
  az_param=""
  if [[ -n "${AWS_AZ}" ]]; then
    az_param="--placement AvailabilityZone=${AWS_AZ}"
    echo "Using specified availability zone: ${AWS_AZ}"
  fi

  SECURITY_GROUP_ID=$(aws ec2 describe-security-groups \
    --region ${AWS_REGION} \
    --filters "Name=group-name,Values=${EC2_SECURITY_GROUP}" \
    --query 'SecurityGroups[0].GroupId' \
    --output text)

  launch_output=$(
    aws ec2 run-instances --region ${AWS_REGION} \
    --image-id ${AMI_IMAGE} \
    --count ${NUM_INSTANCES} \
    --instance-type ${EC2_TYPE} \
    --key-name ${EC2_KEY_NAME} \
    --security-group-ids ${sg_id} \
    --block-device-mappings "${BLK_DEV_JSON}" \
    ${userdata} \
    ${az_param} \
    --output json
  )

  if [[ -z "$launch_output" ]]; then
    echo "Error: Failed to launch instances. Check your AWS credentials and permissions."
    exit 1
  fi

  # Extract all instance IDs and their AZs
  instance_ids=($(echo "$launch_output" | jq -r '.Instances[].InstanceId'))
  azs=($(echo "$launch_output" | jq -r '.Instances[].Placement.AvailabilityZone'))
  
  if [[ ${#instance_ids[@]} -eq 0 ]]; then
    echo "Error: Failed to create instances. $launch_output"
    exit 1
  fi

  # Tag instances with sequential names
  for i in "${!instance_ids[@]}"; do
    aws ec2 create-tags --resources ${instance_ids[$i]} \
        --tags "Key=Name,Value=${INSTANCE_NAME}-$((i+1))"
    echo "Instance ${instance_ids[$i]} launched in availability zone: ${azs[$i]}"
  done
  
  # Set INSTANCE_INFO to first instance ID for compatibility
  INSTANCE_INFO="${instance_ids[0]}"

}

function add_disks() {
  for instance_id in "${instance_ids[@]}"; do
    echo "Adding ${EBS_QTY} ${EBS_TYPE} volumes (${EBS_SIZE}GB) to ${instance_id}"
    ./ebs-create-attach.sh "${instance_id}" "${EBS_TYPE}" "${EBS_SIZE}" "${EBS_QTY}"
  done
}

function register_dns() {
  hosted_zone_id=$(aws route53 list-hosted-zones --query 'HostedZones[0].Id' --output text)
  if [[ -n "$hosted_zone_id" && "$hosted_zone_id" != "None" ]]; then
    for i in "${!public_ips[@]}"; do
      hostname="${INSTANCE_NAME}-$((i+1))"
      fqdn="${hostname}.${DOMAIN}"
      
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
            "Value": "${public_ips[$i]}"
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
        echo "DNS record created: ${fqdn}"
      else
        echo "Failed to create DNS record for ${fqdn}"
      fi
    done
  else
    echo "No Route 53 hosted zone found. Using public IP addresses."
  fi
}

function configure_ceph_nodes() {
    local num_new_instances=${#instance_ids[@]}
    if [[ $num_new_instances -eq 0 ]]; then
        echo "No new instances available to configure for Ceph."
        return
    fi

    local key_source_node_ip=""
    local key_source_node_fqdn_display="" # For logging
    local temp_ceph_pub_key_local="/tmp/ceph_cluster_key_$(date +%s%N).pub" # Unique temp file on control machine

    if [[ -z "$EXISTING_CEPH_ADMIN_FQDN" ]]; then
        # === Initial Cluster Creation Mode ===
        echo "Mode: Initial Ceph Cluster Creation."
        local first_new_node_public_ip="${public_ips[0]}"
        local first_new_node_instance_id="${instance_ids[0]}"
        key_source_node_ip="$first_new_node_public_ip"
        key_source_node_fqdn_display="${INSTANCE_NAME}-1.${DOMAIN}" # Assumes DNS is set or for display

        echo "Bootstrapping Ceph on the first new node: $key_source_node_fqdn_display ($key_source_node_ip)"

        echo "Copying ./ceph-bootstrap.sh to ${EC2_USER}@${key_source_node_ip}:/tmp/ceph-bootstrap.sh..."
        scp -i "${EC2_KEY_FILE}" -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null \
            ./ceph-bootstrap.sh "${EC2_USER}@${key_source_node_ip}:/tmp/ceph-bootstrap.sh"
        if [[ $? -ne 0 ]]; then
            echo "Error: Failed to copy ceph-bootstrap.sh to $key_source_node_fqdn_display. Aborting Ceph configuration."
            return 1
        fi

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
        # === Cluster Expansion Mode ===
        echo "Mode: Expanding existing Ceph Cluster. Admin node for key source: $EXISTING_CEPH_ADMIN_FQDN"
        key_source_node_ip=$(dig +short A "$EXISTING_CEPH_ADMIN_FQDN" | head -n1)
        if [[ -z "$key_source_node_ip" ]]; then
            echo "Error: Could not resolve IP for existing Ceph admin node $EXISTING_CEPH_ADMIN_FQDN. Aborting."
            return 1
        fi
        key_source_node_fqdn_display="$EXISTING_CEPH_ADMIN_FQDN"
        echo "Using existing Ceph admin node $key_source_node_fqdn_display ($key_source_node_ip) as source for Ceph public key."
    fi

    # === Common Steps for Both Modes: Fetch Ceph Public Key ===
    echo "Fetching Ceph public key from $key_source_node_fqdn_display ($key_source_node_ip)..."

    # Prepare Ceph public key on the key_source_node (copy to /tmp, chown, chmod)
    # This requires SSH access to key_source_node_ip using EC2_KEY_FILE (or appropriate auth for existing admin)
    echo "Preparing Ceph public key on $key_source_node_fqdn_display..."
    local prep_key_cmd="sudo cp /etc/ceph/ceph.pub /tmp/ceph.pub && sudo chown ${EC2_USER}:${EC2_USER} /tmp/ceph.pub && sudo chmod 644 /tmp/ceph.pub"
    ssh -i "${EC2_KEY_FILE}" -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null \
        "${EC2_USER}@${key_source_node_ip}" "${prep_key_cmd}"
    if [[ $? -ne 0 ]]; then
        echo "Error: Failed to prepare Ceph public key on $key_source_node_fqdn_display (/tmp/ceph.pub). Aborting key distribution."
        if [[ -z "$EXISTING_CEPH_ADMIN_FQDN" ]]; then # If initial mode, cleanup bootstrap script
            ssh -i "${EC2_KEY_FILE}" -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null \
                "${EC2_USER}@${public_ips[0]}" "rm -f /tmp/ceph-bootstrap.sh"
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
        if [[ -z "$EXISTING_CEPH_ADMIN_FQDN" ]]; then
            ssh -i "${EC2_KEY_FILE}" -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null \
                "${EC2_USER}@${public_ips[0]}" "rm -f /tmp/ceph-bootstrap.sh"
        fi
        return 1
    fi

    # Clean up /tmp/ceph.pub on the key_source_node as it's now on the control machine
    ssh -i "${EC2_KEY_FILE}" -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null \
        "${EC2_USER}@${key_source_node_ip}" "rm -f /tmp/ceph.pub"

    # Determine target nodes for key distribution
    local target_node_indices=()
    if [[ -z "$EXISTING_CEPH_ADMIN_FQDN" ]]; then
        # Initial mode: distribute to new nodes 1 to N-1 (0-indexed: public_ips[1] to public_ips[num_new_instances-1])
        # If only 1 new instance, this loop won't run.
        for i in $(seq 1 $((num_new_instances - 1))); do
            target_node_indices+=($i)
        done
    else
        # Expansion mode: distribute to all new nodes (0 to N-1)
        for i in $(seq 0 $((num_new_instances - 1))); do
            target_node_indices+=($i)
        done
    fi

    if [[ ${#target_node_indices[@]} -gt 0 ]]; then
        echo "Distributing Ceph public key to designated newly created nodes..."
        for i in "${target_node_indices[@]}"; do
            local target_public_ip="${public_ips[$i]}"
            # Instance names are 1-based, array indices are 0-based
            local target_fqdn_display="${INSTANCE_NAME}-$(($i+1)).${DOMAIN}"

            echo "Distributing Ceph public key to $target_fqdn_display ($target_public_ip)..."
            
            # ssh-copy-id from control machine directly to target_public_ip
            # EC2_KEY_FILE is used for authentication.
            local ssh_options_for_copy_id=(
                "-o" "StrictHostKeyChecking=no"
                "-o" "UserKnownHostsFile=/dev/null"
                "-o" "IdentityFile=${EC2_KEY_FILE}"
                "-o" "ConnectTimeout=10"
            )

            if ssh-copy-id -f -i "${temp_ceph_pub_key_local}" "${ssh_options_for_copy_id[@]}" "${EC2_USER}@${target_public_ip}"; then
                echo "Successfully copied Ceph public key to $target_fqdn_display ($target_public_ip)."
            else
                echo "Warning: Failed to copy Ceph public key to $target_fqdn_display ($target_public_ip). Exit status: $?"
            fi
        done
    else
        if [[ -z "$EXISTING_CEPH_ADMIN_FQDN" && $num_new_instances -eq 1 ]]; then
            echo "Only one new instance created in initial mode. Ceph bootstrap done. No further key distribution needed for other *new* nodes."
        elif [[ -n "$EXISTING_CEPH_ADMIN_FQDN" && $num_new_instances -gt 0 ]]; then
             # This case should have been caught by target_node_indices having elements.
             # However, if somehow it's reached and target_node_indices is empty but we are in expansion mode with new instances, it's an anomaly.
             echo "Warning: In expansion mode but no target nodes identified for key distribution. This should not happen if new instances were created."
        else
            echo "No other new nodes to distribute Ceph public key to."
        fi
    fi
    
    # Final Cleanup
    rm -f "${temp_ceph_pub_key_local}" # Remove temp key from control machine
    if [[ -z "$EXISTING_CEPH_ADMIN_FQDN" ]]; then # Initial mode: cleanup bootstrap script from first new node
        ssh -i "${EC2_KEY_FILE}" -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null \
            "${EC2_USER}@${public_ips[0]}" "rm -f /tmp/ceph-bootstrap.sh"
    fi
    echo "Ceph node configuration process complete."
}

function wait_for_instance() {
  public_ips=()  # Reset the global array
  echo "Waiting for instances to be ready..."
  local max_attempts=30
  local attempt=0
  
  while [ $attempt -lt $max_attempts ]; do
      all_ready=true
      for instance_id in "${instance_ids[@]}"; do
          instance_info=$(aws ec2 describe-instances --region ${AWS_REGION} --instance-ids $instance_id \
              --query 'Reservations[0].Instances[0].{State:State.Name,PublicIpAddress:PublicIpAddress,InstanceType:InstanceType}' \
              --output json 2>/dev/null)
          
          if [[ $? -eq 0 ]]; then
              state=$(echo $instance_info | jq -r '.State')
              public_ip=$(echo $instance_info | jq -r '.PublicIpAddress')
              instance_type=$(echo $instance_info | jq -r '.InstanceType')
              
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

  # Collect all public IPs
  for instance_id in "${instance_ids[@]}"; do
      instance_info=$(aws ec2 describe-instances --region ${AWS_REGION} --instance-ids $instance_id \
          --query 'Reservations[0].Instances[0].PublicIpAddress' --output text)
      public_ips+=("$instance_info")
  done

  echo "Waiting for SSH to become available on all instances..."
  for public_ip in "${public_ips[@]}"; do
      max_ssh_attempts=30
      ssh_attempt=0
      sshout=$(mktemp -t ec2-XXX)
      while [ $ssh_attempt -lt $max_ssh_attempts ]; do
          if ssh -i "${EC2_KEY_FILE}" -o ConnectTimeout=5 -o BatchMode=yes -o StrictHostKeyChecking=no ${EC2_USER}@${public_ip} exit 2>${sshout}; then
              echo "SSH is available on $public_ip"
              break
          fi
          printf "Debug ($public_ip): " && cat $sshout
          echo "Waiting for SSH... (Attempt $((ssh_attempt+1))/$max_ssh_attempts)"
          sleep 10
          ssh_attempt=$((ssh_attempt+1))
      done

      if [[ $ssh_attempt -eq $max_ssh_attempts ]]; then
          echo "Error: SSH did not become available for $public_ip within the expected time."
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

launch_instance
wait_for_instance
add_disks
register_dns
configure_ceph_nodes

echo -e "\nInstances created on ${EC2_TYPE} with AMI ${AMI_IMAGE}"
echo "SSH commands to connect:"
FILE2=~${EC2_KEY_FILE#"$HOME"}
for i in "${!public_ips[@]}"; do
    hostnum=$((i+1))
    fqdn="${INSTANCE_NAME}-${hostnum}.${DOMAIN}"
    echo "ssh -i '${FILE2}' ${EC2_USER}@${fqdn}"
done

