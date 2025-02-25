#!/bin/bash

# AWS EC2 Instance Launch Script
# This script launches a c7gd.medium EC2 instance with Rocky Linux 9

# Configuration Variables
AWS_REGION=us-west-2
EC2_TYPE="c7gd.medium"
AMI_IMAGE="ami-03be04a3da3a40226"
ROOT_VOLUME_SIZE=16
REGION="us-west-2"
INSTANCE_NAME="ceph-test-1"
SECURITY_GROUP_ID="sg-021caf5c128f60395"
EC2_KEY_NAME="auto-ec2-405644541454-dirkcli" # Adjust this to your key name
DOMAIN="ai.oregonstate.edu"  # Your domain for DNS record
CLOUD_INIT_FILE="/home/dp/gh/dptests/ec2-cloud-init.txt"  # Path to your cloud-init file
EC2_USER="rocky"

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

  launch_output=$(
    aws ec2 run-instances --region ${AWS_REGION} \
    --image-id ${AMI_IMAGE} \
    --count 1 \
    --instance-type ${EC2_TYPE} \
    --key-name ${EC2_KEY_NAME} \
    --security-group-ids ${sg_id} \
    --block-device-mappings "${BLK_DEV_JSON}" \
    ${userdata} \
    ${az_param} \
    --tag-specifications "ResourceType=instance,Tags=[{Key=Name,Value=${NAME}}]" \
    --query 'Instances[0].{InstanceId:InstanceId,AvailabilityZone:Placement.AvailabilityZone}' \
    --output json
  )

  if [[ -z "$launch_output" ]]; then
    echo "Error: Failed to launch instance. Check your AWS credentials and permissions."
    exit 1
  fi

  # Extract instance ID and AZ from the JSON response
  instance_id=$(echo "$launch_output" | jq -r '.InstanceId')
  az=$(echo "$launch_output" | jq -r '.AvailabilityZone')
  
  if [[ -z "$instance_id" || "$instance_id" == "null" ]]; then
    echo "Error: Failed to create instance. $launch_output"
    exit 1
  fi

  echo "Instance launched in availability zone: ${az}"
  
  # Set INSTANCE_INFO to just the instance ID for compatibility with the rest of the script
  INSTANCE_INFO="$instance_id"
}

function register_dns() {
  # Check for Route 53 hosted zone and register DNS if available
  hosted_zone_id=$(aws route53 list-hosted-zones --query 'HostedZones[0].Id' --output text)
  if [[ -n "$hosted_zone_id" && "$hosted_zone_id" != "None" ]]; then
    domain_name=$(aws route53 get-hosted-zone --id $hosted_zone_id --query 'HostedZone.Name' --output text)
    domain_name=${domain_name%%.} # Remove trailing dot
    fqdn="${NAME}.${domain_name}"
    
    # Store hosted_zone_id and fqdn as global variables
    export HOSTED_ZONE_ID=$hosted_zone_id
    export FQDN=$fqdn
    
    # Create Route 53 DNS record
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

    if aws route53 change-resource-record-sets --hosted-zone-id $hosted_zone_id --change-batch "$change_batch" > /dev/null; then
      echo "DNS record created: ${fqdn}"
      public_dns_name=$fqdn
    else
      echo "Failed to create DNS record. Using public DNS name."
    fi
  else
    echo "No Route 53 hosted zone found. Using public DNS name."
  fi
}

function wait_for_instance() {
  echo "Waiting for instance to be ready..."
  local max_attempts=30
  local attempt=0
  while [ $attempt -lt $max_attempts ]; do
      instance_info=$(aws ec2 describe-instances --region ${AWS_REGION} --instance-ids $instance_id \
          --query 'Reservations[0].Instances[0].{State:State.Name,PublicIpAddress:PublicIpAddress,InstanceType:InstanceType}' \
          --output json 2>/dev/null)
      
      if [[ $? -eq 0 ]]; then
          state=$(echo $instance_info | jq -r '.State')
          public_ip=$(echo $instance_info | jq -r '.PublicIpAddress')
          instance_type=$(echo $instance_info | jq -r '.InstanceType')
          
          if [[ "$state" = "running" && "$public_ip" != "null" && -n "$public_ip" ]]; then
              echo "Instance is ready. Public IP: $public_ip"
              break
          fi
      fi
      
      echo "Waiting for instance to be ready... (Attempt $((attempt+1))/$max_attempts)"
      sleep 10
      attempt=$((attempt+1))
  done

  if [[ $attempt -eq $max_attempts ]]; then
      echo "Error: Instance did not become ready within the expected time."
      exit 1
  fi

  echo "Waiting for SSH to become available..."
  max_ssh_attempts=30
  ssh_attempt=0
  sshout=$(mktemp -t ec2-XXX)
  while [ $ssh_attempt -lt $max_ssh_attempts ]; do
      if ssh -i "${EC2_KEY_FILE}" -o ConnectTimeout=5 -o BatchMode=yes -o StrictHostKeyChecking=no ${EC2_USER}@${public_ip} exit 2>${sshout}; then
          echo "SSH is now available."
          return 0
      fi
      printf "Debug: " && cat $sshout
      echo "Waiting for SSH... (Attempt $((ssh_attempt+1))/$max_ssh_attempts)"
      sleep 10
      ssh_attempt=$((ssh_attempt+1))
  done

  echo "Error: SSH did not become available within the expected time."
  exit 1
}
launch_instance
wait_for_instance
register_dns

echo "Starting $DNS_NAME ($INSTANCE_ID) on $INSTANCE_TYPE with $AMI_NAME"
echo "Run one of these SSH commands to connect:"
echo "ssh -i ~/.ssh/${KEY_NAME}.pem ${EC2_USER}@${DNS_NAME}"
echo "ssh -i ~/.ssh/${KEY_NAME}.pem ${EC2_USER}@${PUBLIC_IP}"
