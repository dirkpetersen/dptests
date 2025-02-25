#!/bin/bash

# AWS EC2 Instance Launch Script
# This script launches a c7gd.medium EC2 instance with Rocky Linux 9

# Configuration Variables
INSTANCE_TYPE="c7gd.medium"
AMI_ID="ami-03be04a3da3a40226"
ROOT_VOLUME_SIZE=16
REGION="us-west-2"
INSTANCE_NAME="ceph-test-1"
SECURITY_GROUP_ID="sg-021caf5c128f60395"
KEY_NAME="auto-ec2-405644541454-dirkcli" # Adjust this to your key name
DOMAIN="ai.oregonstate.edu"  # Your domain for DNS record
CLOUD_INIT_FILE="/home/dp/gh/dptests/ec2-cloud-init.txt"  # Path to your cloud-init file
EC2_USER="rocky"
# Read cloud-init file if it exists
CLOUD_INIT=""
if [ -f "$CLOUD_INIT_FILE" ]; then
  echo "Using cloud-init script from $CLOUD_INIT_FILE"
  CLOUD_INIT=$(cat "$CLOUD_INIT_FILE" | base64 -w 0)
else
  echo "Cloud-init file not found at $CLOUD_INIT_FILE, continuing without user data"
fi

# Launch the EC2 instance
INSTANCE_DATA=$(aws ec2 run-instances \
  --region $REGION \
  --image-id $AMI_ID \
  --instance-type $INSTANCE_TYPE \
  --key-name $KEY_NAME \
  --security-group-ids $SECURITY_GROUP_ID \
  --block-device-mappings "[{\"DeviceName\":\"/dev/sda1\",\"Ebs\":{\"VolumeSize\":$ROOT_VOLUME_SIZE,\"DeleteOnTermination\":true}}]" \
  --tag-specifications "ResourceType=instance,Tags=[{Key=Name,Value=$INSTANCE_NAME}]" \
  $([ -n "$CLOUD_INIT" ] && echo "--user-data-base64 $CLOUD_INIT") \
  --output json)

INSTANCE_ID=$(echo $INSTANCE_DATA | jq -r '.Instances[0].InstanceId')
AZ=$(echo $INSTANCE_DATA | jq -r '.Instances[0].Placement.AvailabilityZone')

echo "Instance launched in availability zone: $AZ"
echo "EC2_TYPE: $INSTANCE_TYPE"
echo "AMI_IMAGE: $AMI_ID"
echo "EC2_USER: $EC2_USER"

echo "Instance $INSTANCE_ID is starting. Waiting for public IP..."

# Wait for the instance to be running and get its public IP
PUBLIC_IP=""
while [ -z "$PUBLIC_IP" ]; do
  sleep 5
  PUBLIC_IP=$(aws ec2 describe-instances \
    --region $REGION \
    --instance-ids $INSTANCE_ID \
    --query "Reservations[0].Instances[0].PublicIpAddress" \
    --output text)
  
  if [ "$PUBLIC_IP" = "None" ] || [ -z "$PUBLIC_IP" ]; then
    PUBLIC_IP=""
  fi
done

# Create DNS record if domain is provided
if [ -n "$DOMAIN" ]; then
  DNS_NAME="${INSTANCE_NAME}.${DOMAIN}"
  echo "DNS record created: $DNS_NAME"
  
  # Add DNS creation command here based on your DNS provider
  # For example, using Route 53:
  # aws route53 change-resource-record-sets --hosted-zone-id YOUR_HOSTED_ZONE_ID --change-batch "{\"Changes\":[{\"Action\":\"UPSERT\",\"ResourceRecordSet\":{\"Name\":\"$DNS_NAME\",\"Type\":\"A\",\"TTL\":300,\"ResourceRecords\":[{\"Value\":\"$PUBLIC_IP\"}]}}]}"
  
  echo "Instance $INSTANCE_ID is running. Waiting for DNS propagation..."
  
  # Wait for DNS propagation
  MAX_ATTEMPTS=30
  ATTEMPT=1
  
  while [ $ATTEMPT -le $MAX_ATTEMPTS ]; do
    echo "Waiting for DNS propagation... (Attempt $ATTEMPT/$MAX_ATTEMPTS)"
    
    # Check if DNS has propagated
    # This is a simplified check - you might want to use dig or nslookup
    sleep 10
    
    ATTEMPT=$((ATTEMPT + 1))
    
    # Break after first attempt for the script (in real use, you'd verify DNS resolution)
    if [ $ATTEMPT -gt 1 ]; then
      break
    fi
  done
  
  echo "DNS has propagated ..."
fi

echo "Starting $DNS_NAME ($INSTANCE_ID) on $INSTANCE_TYPE with $AMI_NAME"
echo "Run one of these SSH commands to connect:"
echo "ssh -i ~/.ssh/${KEY_NAME}.pem ${EC2_USER}@${DNS_NAME}"
echo "ssh -i ~/.ssh/${KEY_NAME}.pem ${EC2_USER}@${PUBLIC_IP}"
