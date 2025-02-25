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
DEVICE_LETTERS=(b c d e f g h i j k l m n o p q r s t u v w x y z)
COUNTER=2

# Detect whether the instance uses NVMe or Xen
IS_NITRO=$(aws ec2 describe-instances --instance-ids "$INSTANCE_ID" --query "Reservations[0].Instances[0].EnaSupport" --output text 2>/dev/null)

for VOLUME_ID in "${VOLUME_IDS[@]}"; do
    DEVICE_NAME="/dev/xvd${DEVICE_LETTERS[$COUNTER]}"
    if [[ "$IS_NITRO" != "True" ]]; then
        # Xen instances use xvdX naming
         echo "Attaching volume $VOLUME_ID to instance $INSTANCE_ID as $DEVICE_NAME"
    fi
    
    aws ec2 attach-volume \
        --instance-id "$INSTANCE_ID" \
        --volume-id "$VOLUME_ID" \
        --device "$DEVICE_NAME"

    ((COUNTER++))
done

echo "All volumes created and attached successfully."

