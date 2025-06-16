#! /bin/bash

CEPH_RELEASE=19.2.2  # replace this with the active release
#CEPH_RELEASE=18.2.7

# Configuration: Number of HDDs that should share one SSD for DB
: "${HDDS_PER_SSD:=6}"

# Accept MON IP as first argument (for joining existing cluster)
MON_IP_ARG="$1"

curl -o /usr/local/bin/cephadm https://download.ceph.com/rpm-${CEPH_RELEASE}/el9/noarch/cephadm
chmod +x /usr/local/bin/cephadm
# Try to get IMDSv2 token with a short timeout to see if we are on EC2
TOKEN=$(curl --connect-timeout 2 -X PUT "http://169.254.169.254/latest/api/token" -H "X-aws-ec2-metadata-token-ttl-seconds: 21600" 2>/dev/null)
# Check if token was successfully retrieved
if [[ -n "$TOKEN" ]]; then
    echo "Running on EC2, using metadata service version 2..."
    # Get internal IP using the token
    export INTERNAL_IP=$(curl -H "X-aws-ec2-metadata-token: $TOKEN" -s http://169.254.169.254/latest/meta-data/local-ipv4)
else
    # If the token retrieval fails, we will fall back to using the hostname command
    echo "EC2 metadata service not available, using hostname command..."
    export INTERNAL_IP=$(hostname -I | awk '{print $1}')
fi

# Configure udev rules for AWS EBS volumes before deploying Ceph
echo "Configuring udev rules for AWS EBS volumes..."
cat > /etc/udev/rules.d/99-aws-ebs.rules << 'EOF'
# Add this content (targets ONLY "Amazon Elastic Block Store"):
SUBSYSTEM=="block", KERNEL=="nvme*", ENV{ID_MODEL}=="Amazon Elastic Block Store", ENV{DEVTYPE}=="disk", ATTR{queue/rotational}="1"
EOF

# Reload udev rules
echo "Reloading udev rules..."
udevadm control --reload-rules
udevadm trigger

# Determine mode: bootstrap new cluster or join existing cluster
if [ -n "$MON_IP_ARG" ]; then
    echo "Using provided MON IP: $MON_IP_ARG (joining existing cluster)"
    MON_IP="$MON_IP_ARG"
    
    # For joining existing cluster, we need to get the cluster configuration first
    echo "Joining existing Ceph cluster at $MON_IP..."
    
    # Try to get cluster config from MON node
    echo "Attempting to retrieve cluster configuration from MON at $MON_IP..."
    # For now, we'll just install cephadm and set up for OSD creation
    # The actual cluster joining will be handled by the orchestrator from the MON node
    
    echo "Ceph cluster join mode - OSDs will be added by orchestrator from MON node"
    # Skip bootstrap, just prepare for OSD creation
    SKIP_BOOTSTRAP=true
else
    echo "No MON IP provided, using internal IP: $INTERNAL_IP (bootstrapping new cluster)"
    MON_IP="$INTERNAL_IP"
    SKIP_BOOTSTRAP=false
fi

if [ "$SKIP_BOOTSTRAP" != "true" ]; then
    #/usr/local/bin/cephadm bootstrap --allow-fqdn-hostname --mon-ip ${MON_IP} --ssh-user root
    /usr/local/bin/cephadm bootstrap --mon-ip ${MON_IP} --ssh-user root

    # Restart ceph services after bootstrap
    echo "Restarting Ceph services..."
    systemctl restart ceph.target

    # Wait for Ceph cluster to be ready
    echo "Waiting for Ceph cluster to be ready..."
    sleep 30
else
    echo "Skipping bootstrap - this node will join the existing cluster via orchestrator"
    # Just ensure cephadm is available for when the orchestrator adds this node
    exit 0
fi

# # Wait for bootstrap keyring to be available
# echo "Waiting for bootstrap keyring..."
# while [ ! -f /var/lib/ceph/bootstrap-osd/ceph.keyring ]; do
#     echo "Bootstrap keyring not found, waiting..."
#     sleep 10
# done
# echo "Bootstrap keyring found, proceeding with OSD creation..."

# Create OSDs with shared WAL and DB on SSD
echo "Creating OSDs with shared WAL/DB configuration..."

# Use ceph orch device ls to detect available devices  
echo "Detecting available storage devices using ceph orch device ls..."
SSD_DEVICES=()
HDD_DEVICES=()

# Get available devices from ceph orch device ls
device_output=$(/usr/local/bin/cephadm shell --fsid $FSID -c /etc/ceph/ceph.conf -k /etc/ceph/ceph.client.admin.keyring -- \
    ceph orch device ls --format json 2>/dev/null)

if [ $? -eq 0 ] && [ -n "$device_output" ]; then
    echo "Using ceph orch device ls to detect devices..."
    # Parse JSON output to get available devices
    available_devices=$(echo "$device_output" | python3 -c "
import json, sys
try:
    data = json.load(sys.stdin)
    for host_data in data:
        for device in host_data.get('devices', []):
            if device.get('available', False) and not device.get('rejected_reasons', []):
                path = device.get('path', '')
                rotational = device.get('rotational', True) 
                if path:
                    print(f'{path}:{rotational}')
except Exception as e:
    print(f'Error parsing device list: {e}', file=sys.stderr)
")
    
    # Process the available devices
    while IFS=':' read -r device rotational; do
        if [ -n "$device" ]; then
            if [ "$rotational" = "False" ]; then
                SSD_DEVICES+=("$device")
                echo "Found available SSD: $device"
            else
                HDD_DEVICES+=("$device") 
                echo "Found available HDD: $device"
            fi
        fi
    done <<< "$available_devices"
else
    echo "ceph orch device ls failed, falling back to manual detection..."
    # Fallback to manual detection for block devices only (not partitions)
    for device in /dev/sd* /dev/nvme*n*; do
        if [ -b "$device" ] && [[ ! "$device" =~ p[0-9]+$ ]]; then
            # Check if device is truly available (no filesystem, not mounted, not in use)
            if ! lsblk -n -o FSTYPE "$device" 2>/dev/null | grep -q .; then
                if ! findmnt -S "$device" >/dev/null 2>&1; then
                    # Check rotational attribute (0=SSD, 1=HDD)
                    rotational=$(cat /sys/block/$(basename $device)/queue/rotational 2>/dev/null || echo "1")
                    if [ "$rotational" = "0" ]; then
                        SSD_DEVICES+=("$device")
                        echo "Found available SSD: $device"
                    else
                        HDD_DEVICES+=("$device")
                        echo "Found available HDD: $device"
                    fi
                fi
            fi
        fi
    done
fi

# Check if we have devices to work with
if [ ${#HDD_DEVICES[@]} -eq 0 ] && [ ${#SSD_DEVICES[@]} -eq 0 ]; then
    echo "No available devices found for OSD creation."
    exit 0
fi

# Get cluster info
FSID=$(cat /etc/ceph/ceph.conf | grep fsid | awk '{print $3}')
HOSTNAME=$(hostname)

# Check if we have at least one SSD for WAL/DB
if [ ${#SSD_DEVICES[@]} -eq 0 ]; then
    echo "No SSD devices found. Creating OSDs with data-only configuration..."
    for hdd in "${HDD_DEVICES[@]}"; do
        echo "Creating OSD on $hdd (data only)..."
        /usr/local/bin/cephadm shell --fsid $FSID -c /etc/ceph/ceph.conf -k /etc/ceph/ceph.client.admin.keyring -- \
          ceph orch daemon add osd "$HOSTNAME:data_devices=$hdd"
        if [ $? -ne 0 ]; then
            echo "Warning: Failed to create OSD on $hdd"
        fi
    done
else
    echo "Creating OSDs with shared DB configuration (up to $HDDS_PER_SSD HDDs per SSD)..."
    
    # Process HDDs in groups, assigning one SSD for DB per group
    hdd_index=0
    ssd_index=0
    
    while [ $hdd_index -lt ${#HDD_DEVICES[@]} ] && [ $ssd_index -lt ${#SSD_DEVICES[@]} ]; do
        # Calculate how many HDDs to assign to this SSD (respect HDDS_PER_SSD limit)
        remaining_hdds=$((${#HDD_DEVICES[@]} - hdd_index))
        hdds_for_this_ssd=$((remaining_hdds < HDDS_PER_SSD ? remaining_hdds : HDDS_PER_SSD))
        
        # Build comma-separated list of HDDs for this batch
        data_devices_list=""
        for ((i=0; i<hdds_for_this_ssd; i++)); do
            if [ $((hdd_index + i)) -lt ${#HDD_DEVICES[@]} ]; then
                if [ -n "$data_devices_list" ]; then
                    data_devices_list="$data_devices_list,${HDD_DEVICES[$((hdd_index + i))]}"
                else
                    data_devices_list="${HDD_DEVICES[$((hdd_index + i))]}"
                fi
            fi
        done
        
        if [ -n "$data_devices_list" ]; then
            ssd_device="${SSD_DEVICES[$ssd_index]}"
            echo "Creating $hdds_for_this_ssd OSDs with data on [$data_devices_list], DB on $ssd_device..."
            
            /usr/local/bin/cephadm shell --fsid $FSID -c /etc/ceph/ceph.conf -k /etc/ceph/ceph.client.admin.keyring -- \
              ceph orch daemon add osd "$HOSTNAME:data_devices=$data_devices_list,db_devices=$ssd_device,osds_per_device=1"
            
            if [ $? -ne 0 ]; then
                echo "Warning: Failed to create OSDs with data devices [$data_devices_list] and DB on $ssd_device"
            fi
        fi
        
        # Move to next group
        hdd_index=$((hdd_index + hdds_for_this_ssd))
        ssd_index=$((ssd_index + 1))
    done
    
    # Use any remaining SSDs for data-only OSDs
    while [ $ssd_index -lt ${#SSD_DEVICES[@]} ]; do
        ssd="${SSD_DEVICES[$ssd_index]}"
        echo "Creating SSD OSD on $ssd (data only)..."
        /usr/local/bin/cephadm shell --fsid $FSID -c /etc/ceph/ceph.conf -k /etc/ceph/ceph.client.admin.keyring -- \
          ceph orch daemon add osd "$HOSTNAME:data_devices=$ssd"
        if [ $? -ne 0 ]; then
            echo "Warning: Failed to create SSD OSD on $ssd"
        fi
        ssd_index=$((ssd_index + 1))
    done
fi

echo "OSD creation complete. Cluster status:"
/usr/local/bin/cephadm shell -- ceph -s
