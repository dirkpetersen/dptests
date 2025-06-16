#! /bin/bash

CEPH_RELEASE=19.2.2  # replace this with the active release
#CEPH_RELEASE=18.2.7

# Configuration: Number of HDDs that should share one SSD for DB
: "${HDDS_PER_SSD:=6}"

# No arguments needed - script auto-detects what to do

curl -o /usr/local/bin/cephadm https://download.ceph.com/rpm-${CEPH_RELEASE}/el9/noarch/cephadm
chmod +x /usr/local/bin/cephadm
# Try to get IMDSv2 token with a short timeout to see if we are on EC2
TOKEN=$(curl --connect-timeout 2 -X PUT "http://169.254.169.254/latest/api/token" -H "X-aws-ec2-metadata-token-ttl-seconds: 21600" 2>/dev/null)
# Check if token was successfully retrieved
if [[ -n "${TOKEN}" ]]; then
    echo "Running on EC2, using metadata service version 2..."
    # Get internal IP using the token
    export INTERNAL_IP=$(curl -H "X-aws-ec2-metadata-token: ${TOKEN}" -s http://169.254.169.254/latest/meta-data/local-ipv4)
else
    # If the token retrieval fails, we will fall back to using the hostname command
    echo "EC2 metadata service not available, using hostname command..."
    export INTERNAL_IP=$(hostname -I | awk '{print $1}')
fi

# Wait for podman to be installed by cloud-init
echo "Waiting for podman to be available (installed by cloud-init)..."
podman_check_attempts=30 # Approx 5 minutes (30 * 10s)
podman_check_count=0
while ! command -v podman >/dev/null 2>&1; do
    podman_check_count=$((podman_check_count + 1))
    if [[ ${podman_check_count} -ge ${podman_check_attempts} ]]; then
        echo "Error: podman did not become available after ${podman_check_attempts} attempts."
        exit 1
    fi
    echo "podman not yet available, waiting 10s... (Attempt ${podman_check_count}/${podman_check_attempts})"
    sleep 10
done
echo "podman is available."

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

# Auto-detect what to do based on cluster state
echo "Auto-detecting cluster state..."

# Check if this node already has Ceph installed
if [[ -f /etc/ceph/ceph.conf ]]; then
    echo "Ceph already configured on this node"
    
    # Check if OSDs already exist on this node
    existing_osds=$(ls /var/lib/ceph/osd/ 2>/dev/null | wc -l)
    if [[ ${existing_osds} -gt 0 ]]; then
        echo "Found ${existing_osds} existing OSDs on this node. Nothing to do."
        echo "Ceph cluster status:"
        /usr/local/bin/cephadm shell -- ceph -s 2>/dev/null || echo "Could not get cluster status"
        exit 0
    else
        echo "No OSDs found - will create OSDs only"
        SKIP_BOOTSTRAP=true
        SKIP_OSD_CREATION=false
    fi
else 
    # Check if we can find an existing cluster in our subnet
    echo "Scanning subnet for existing Ceph cluster..."
    
    # Get our subnet (first 3 octets of internal IP)
    SUBNET_PREFIX=$(echo "${INTERNAL_IP}" | cut -d. -f1-3)
    EXISTING_MON_IP=""
    
    # Scan common IP ranges in our subnet for existing Ceph MON
    for last_octet in $(seq 1 254); do
        candidate_ip="${SUBNET_PREFIX}.${last_octet}"
        
        # Skip our own IP
        if [[ "${candidate_ip}" == "${INTERNAL_IP}" ]]; then
            continue
        fi
        
        # Try to connect to Ceph MON port (quick check)
        if timeout 2 nc -z "${candidate_ip}" 6789 2>/dev/null; then
            echo "Found potential Ceph cluster at ${candidate_ip}"
            EXISTING_MON_IP="${candidate_ip}"
            break
        fi
    done
    
    if [[ -n "${EXISTING_MON_IP}" ]]; then
        echo "Joining existing Ceph cluster at ${EXISTING_MON_IP}"
        MON_IP="${EXISTING_MON_IP}"
        SKIP_BOOTSTRAP=true
        SKIP_OSD_CREATION=false
        
        # Try to get cluster config from existing MON
        echo "Attempting to get cluster configuration from ${EXISTING_MON_IP}..."
        # We'll let the orchestrator handle adding this node
    else
        echo "No existing cluster found - bootstrapping new cluster with IP: ${INTERNAL_IP}"
        MON_IP="${INTERNAL_IP}"
        SKIP_BOOTSTRAP=false
        SKIP_OSD_CREATION=false
    fi

if [[ "${SKIP_BOOTSTRAP}" == "false" ]]; then
    #/usr/local/bin/cephadm bootstrap --allow-fqdn-hostname --mon-ip ${MON_IP} --ssh-user root
    /usr/local/bin/cephadm bootstrap --mon-ip "${MON_IP}" --ssh-user root

    # Restart ceph services after bootstrap
    echo "Restarting Ceph services..."
    systemctl restart ceph.target

    # Wait for Ceph cluster to be ready
    echo "Waiting for Ceph cluster to be ready..."
    sleep 30
else
    echo "Skipping bootstrap - this node will join the existing cluster"
    
    if [[ -n "${EXISTING_MON_IP}" ]]; then
        # Try to get cluster config and join
        echo "Attempting to join cluster at ${EXISTING_MON_IP}..."
        echo "Adding this node to the existing cluster via orchestrator..."
        
        # The first node in the cluster should add this node
        # We'll use SSH to execute the orchestrator command on the MON node
        # This requires the Ceph SSH key to be already distributed
        
        HOSTNAME="$(hostname)"
        
        # Check if we can SSH to the MON node to add ourselves
        if [[ -x "$(command -v ssh)" ]]; then
            echo "Attempting to add this node to cluster via SSH to MON node..."
            
            # Try to add this host to the cluster from the MON node
            # This assumes SSH keys are set up (which ec2-create-instances.sh should handle)
            if ssh -o ConnectTimeout=10 -o StrictHostKeyChecking=no "root@${EXISTING_MON_IP}" \
                "/usr/local/bin/cephadm shell -- ceph orch host add ${HOSTNAME} ${INTERNAL_IP}" 2>/dev/null; then
                echo "Successfully added ${HOSTNAME} to cluster."
                # Now we can create OSDs
                SKIP_OSD_CREATION=false
            else
                echo "Could not automatically add to cluster. Manual addition required:"
                echo "Run on MON node: ceph orch host add ${HOSTNAME} ${INTERNAL_IP}"
                SKIP_OSD_CREATION=true
            fi
        else
            echo "SSH not available. Manual cluster joining required:"
            echo "Run on MON node: ceph orch host add ${HOSTNAME} ${INTERNAL_IP}"
            SKIP_OSD_CREATION=true
        fi
    fi
fi

# # Wait for bootstrap keyring to be available
# echo "Waiting for bootstrap keyring..."
# while [ ! -f /var/lib/ceph/bootstrap-osd/ceph.keyring ]; do
#     echo "Bootstrap keyring not found, waiting..."
#     sleep 10
# done
# echo "Bootstrap keyring found, proceeding with OSD creation..."

# Create OSDs with shared WAL and DB on SSD (only if not skipping)
if [[ "${SKIP_OSD_CREATION}" != "true" ]]; then
    echo "Creating OSDs with shared WAL/DB configuration..."
else
    echo "Skipping OSD creation for this mode."
    exit 0
fi

# Use ceph orch device ls to detect available devices  
echo "Detecting available storage devices using ceph orch device ls..."
SSD_DEVICES=()
HDD_DEVICES=()

# Get cluster info
FSID="$(grep fsid /etc/ceph/ceph.conf | awk '{print $3}')"
HOSTNAME="$(hostname)"

# Get available devices from ceph orch device ls
device_output=$(/usr/local/bin/cephadm shell --fsid "${FSID}" -c /etc/ceph/ceph.conf -k /etc/ceph/ceph.client.admin.keyring -- \
    ceph orch device ls --format json 2>/dev/null)

if [[ $? -eq 0 && -n "${device_output}" ]]; then
    echo "Using ceph orch device ls to detect devices..."
    
    available_devices=$(echo "${device_output}" | python3 -c "
import json, sys
try:
    data = json.load(sys.stdin)
    # Handle both possible structures
    if isinstance(data, list):
        # If it's a list of host objects
        for host_data in data:
            if isinstance(host_data, dict) and 'devices' in host_data:
                for device in host_data['devices']:
                    if device.get('available', False) and not device.get('rejected_reasons', []):
                        path = device.get('path', '')
                        # Check for human_readable_type first, then rotational field
                        if 'human_readable_type' in device:
                            # Map human_readable_type to rotational boolean
                            device_type = device.get('human_readable_type', '').lower()
                            rotational = device_type not in ['ssd', 'nvme']
                        elif 'rotational' in device:
                            rotational = device.get('rotational', True)
                        else:
                            rotational = True  # Default to HDD if unknown
                        if path:
                            print(f'{path}:{rotational}')
            # If it's a direct list of devices
            elif isinstance(host_data, dict) and 'path' in host_data:
                if host_data.get('available', False) and not host_data.get('rejected_reasons', []):
                    path = host_data.get('path', '')
                    # Check for human_readable_type first, then rotational field
                    if 'human_readable_type' in host_data:
                        # Map human_readable_type to rotational boolean
                        device_type = host_data.get('human_readable_type', '').lower()
                        rotational = device_type not in ['ssd', 'nvme']
                    elif 'rotational' in host_data:
                        rotational = host_data.get('rotational', True)
                    else:
                        rotational = True  # Default to HDD if unknown
                    if path:
                        print(f'{path}:{rotational}')
except Exception as e:
    print(f'Error parsing device list: {e}', file=sys.stderr)
")
    
    # Process the available devices
    while IFS=':' read -r device rotational; do
        if [[ -n "${device}" ]]; then
            if [[ "${rotational}" == "False" ]]; then
                SSD_DEVICES+=("${device}")
                echo "Found available SSD: ${device}"
            else
                HDD_DEVICES+=("${device}") 
                echo "Found available HDD: ${device}"
            fi
        fi
    done <<< "${available_devices}"
else
    echo "Error: Failed to get device list from ceph orchestrator. Cannot proceed with OSD creation."
    echo "Please ensure the Ceph cluster is healthy and this node is properly joined."
    exit 1
fi

# Check if we have devices to work with
if [[ ${#HDD_DEVICES[@]} -eq 0 && ${#SSD_DEVICES[@]} -eq 0 ]]; then
    echo "No available devices found for OSD creation."
    exit 0
fi

# Check if we have at least one SSD for WAL/DB
if [[ ${#SSD_DEVICES[@]} -eq 0 ]]; then
    echo "No SSD devices found. Creating OSDs with data-only configuration..."
    for hdd in "${HDD_DEVICES[@]}"; do
        echo "Creating OSD on ${hdd} (data only)..."
        /usr/local/bin/cephadm shell --fsid "${FSID}" -c /etc/ceph/ceph.conf -k /etc/ceph/ceph.client.admin.keyring -- \
          ceph orch daemon add osd "${HOSTNAME}:data_devices=${hdd}"
        if [[ $? -ne 0 ]]; then
            echo "Warning: Failed to create OSD on ${hdd}"
        fi
    done
else
    echo "Creating OSDs with shared DB configuration (up to ${HDDS_PER_SSD} HDDs per SSD)..."
    
    # Process HDDs in groups, assigning one SSD for DB per group
    hdd_index=0
    ssd_index=0
    
    while [[ ${hdd_index} -lt ${#HDD_DEVICES[@]} && ${ssd_index} -lt ${#SSD_DEVICES[@]} ]]; do
        # Calculate how many HDDs to assign to this SSD (respect HDDS_PER_SSD limit)
        remaining_hdds=$((${#HDD_DEVICES[@]} - hdd_index))
        hdds_for_this_ssd=$((remaining_hdds < HDDS_PER_SSD ? remaining_hdds : HDDS_PER_SSD))
        
        # Build comma-separated list of HDDs for this batch
        data_devices_list=""
        for ((i=0; i<hdds_for_this_ssd; i++)); do
            if [[ $((hdd_index + i)) -lt ${#HDD_DEVICES[@]} ]]; then
                if [[ -n "${data_devices_list}" ]]; then
                    data_devices_list="$data_devices_list,${HDD_DEVICES[$((hdd_index + i))]}"
                else
                    data_devices_list="${HDD_DEVICES[$((hdd_index + i))]}"
                fi
            fi
        done
        
        if [[ -n "${data_devices_list}" ]]; then
            ssd_device="${SSD_DEVICES[$ssd_index]}"
            echo "Creating $hdds_for_this_ssd OSDs with data on [$data_devices_list], DB on $ssd_device..."
            
            /usr/local/bin/cephadm shell --fsid "${FSID}" -c /etc/ceph/ceph.conf -k /etc/ceph/ceph.client.admin.keyring -- \
              ceph orch daemon add osd "$HOSTNAME:data_devices=$data_devices_list,db_devices=$ssd_device,osds_per_device=1"
            
            if [[ $? -ne 0 ]]; then
                echo "Warning: Failed to create OSDs with data devices [$data_devices_list] and DB on $ssd_device"
            fi
        fi
        
        # Move to next group
        hdd_index=$((hdd_index + hdds_for_this_ssd))
        ssd_index=$((ssd_index + 1))
    done
    
    # Use any remaining SSDs for data-only OSDs
    while [[ ${ssd_index} -lt ${#SSD_DEVICES[@]} ]]; do
        ssd="${SSD_DEVICES[$ssd_index]}"
        echo "Creating SSD OSD on $ssd (data only)..."
        /usr/local/bin/cephadm shell --fsid $FSID -c /etc/ceph/ceph.conf -k /etc/ceph/ceph.client.admin.keyring -- \
          ceph orch daemon add osd "$HOSTNAME:data_devices=$ssd"
        if [[ $? -ne 0 ]]; then
            echo "Warning: Failed to create SSD OSD on $ssd"
        fi
        ssd_index=$((ssd_index + 1))
    done
fi

echo "OSD creation complete on ${HOSTNAME}. Cluster status:"
/usr/local/bin/cephadm shell -- ceph -s
fi
