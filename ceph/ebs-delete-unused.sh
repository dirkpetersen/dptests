#! /bin/bash

for VOL in $(aws ec2 describe-volumes --filters "Name=status,Values=available" --query "Volumes[*].VolumeId" --output text); do     
  echo "Deleting volume: $VOL"     
  aws ec2 delete-volume --volume-id "$VOL"
done

