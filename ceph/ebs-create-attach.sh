#!/bin/bash

# Ensure all required arguments are provided
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <instance_id> <volume_type> <size> <quantity>"
    exit 1
fi

INSTANCE_ID=$1
VOLUME_TYPE=$2
SIZE=$3
QUANTITY=$4

# Get the Availability Zone of the instance
AVAILABILITY_ZONE=$(aws ec2 describe-instances --instance-ids "$INSTANCE_ID" --query "Reservations[0].Instances[0].Placement.AvailabilityZone" --output text)

# Ensure instance exists
if [ -z "$AVAILABILITY_ZONE" ]; then
    echo "Error: Instance $INSTANCE_ID not found."
    exit 1
fi

echo "Creating $QUANTITY EBS volumes of type $VOLUME_TYPE, size ${SIZE}GB in $AVAILABILITY_ZONE"

VOLUME_IDS=()

# Create EBS Volumes
for i in $(seq 1 $QUANTITY); do
    VOLUME_ID=$(aws ec2 create-volume \
        --availability-zone "$AVAILABILITY_ZONE" \
        --volume-type "$VOLUME_TYPE" \
        --size "$SIZE" \
        --tag-specifications "ResourceType=volume,Tags=[{Key=Name,Value=ExtraVolume$i}]" \
        --query "VolumeId" --output text)

    if [ -z "$VOLUME_ID" ]; then
        echo "Error: Failed to create volume $i."
        exit 1
    fi

    echo "Created volume $VOLUME_ID"
    VOLUME_IDS+=("$VOLUME_ID")
done

# Attach volumes to the instance
# For Nitro-based instances, EBS volumes are exposed to the OS as NVMe devices (e.g., /dev/nvme0n1, /dev/nvme1n1).
# The --device parameter for the 'aws ec2 attach-volume' command acts as a hint.
# AWS recommends using hints like /dev/sd[f-p] for Nitro instances and /dev/xvd[f-p] for Xen instances.
# This script uses these recommended hints. The OS on a Nitro instance will still see the devices as /dev/nvmeXnY.
DEVICE_LETTERS=(f g h i j k l m n o p q r s t u v w x y z) # Start from 'f' to avoid common root/ephemeral device letters
COUNTER=0 # Initialize counter for the DEVICE_LETTERS array

# Detect whether the instance uses NVMe (Nitro) based on EnaSupport, a reliable indicator.
IS_NITRO=$(aws ec2 describe-instances --instance-ids "$INSTANCE_ID" --query "Reservations[0].Instances[0].EnaSupport" --output text 2>/dev/null)

for VOLUME_ID in "${VOLUME_IDS[@]}"; do
    device_hint_prefix=""
    if [[ "$IS_NITRO" == "True" ]]; then
        device_hint_prefix="/dev/sd"
    else
        device_hint_prefix="/dev/xvd"
    fi
    actual_device_hint="${device_hint_prefix}${DEVICE_LETTERS[$COUNTER]}"
    
    echo "Attaching volume $VOLUME_ID to instance $INSTANCE_ID with device hint $actual_device_hint"
    
    aws ec2 attach-volume \
        --instance-id "$INSTANCE_ID" \
        --volume-id "$VOLUME_ID" \
        --device "$actual_device_hint"

    ((COUNTER++))
done

echo "All volumes created and attached successfully."

