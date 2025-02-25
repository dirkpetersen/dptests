#!/bin/bash

declare -a public_ips

# AWS EC2 Instance Launch Script
# This script launches a c7gd.medium EC2 instance with Rocky Linux 9

# Get number of instances from command line
NUM_INSTANCES=${1:-1}
if [[ ! "$NUM_INSTANCES" =~ ^[1-9][0-9]*$ ]]; then
    echo "Error: Number of instances must be a positive integer"
    exit 1
fi

# Configuration Variables
AWS_REGION=us-west-2
EC2_TYPE="c7gd.medium"
AMI_IMAGE="ami-03be04a3da3a40226"
ROOT_VOLUME_SIZE=16
INSTANCE_NAME="ceph-test"
DOMAIN="ai.oregonstate.edu"  # Your domain for DNS record
CLOUD_INIT_FILE="ec2-cloud-init.txt"  # Path to your cloud-init file
EC2_USER="rocky"
: "${EC2_SECURITY_GROUP:="SSH-HTTP-ICMP"}"

function launch_instance() {
  # Check for ec2-cloud-init.txt file
  userdata=""
  if [[ -n "$CLOUD_INIT_FILE" ]]; then  
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

prepare_ssh_keys
launch_instance
wait_for_instance
register_dns

function cleanup() {
  rm -f "$ROOT_SSH_KEY_FILE" "${ROOT_SSH_KEY_FILE}.pub"
}
trap cleanup EXIT

echo -e "\nInstances created on ${EC2_TYPE} with AMI ${AMI_IMAGE}"
echo "SSH commands to connect:"
FILE2=~${EC2_KEY_FILE#"$HOME"}
for i in "${!public_ips[@]}"; do
    hostnum=$((i+1))
    fqdn="${INSTANCE_NAME}-${hostnum}.${DOMAIN}"
    echo "ssh -i '${FILE2}' ${EC2_USER}@${fqdn}"
done

