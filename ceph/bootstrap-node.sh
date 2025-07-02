#! /bin/bash

CEPH_RELEASE=19.2.2  # replace this with the active release
#CEPH_RELEASE=18.2.7

# Configuration: Number of HDDs that should share one SSD for DB
: "${HDDS_PER_SSD:=2}"

# Accept mode as first argument: 'first' or 'others'
MODE="$1"

if [[ "$MODE" != "first" && "$MODE" != "others" ]]; then
    echo "Usage: $0 {first|others}"
    echo "  first  - Bootstrap the first node (MON node) of a new Ceph cluster"
    echo "  others - Join an existing Ceph cluster as an additional node"
    exit 1
fi

# Function to wait for cluster integration on nodes joining existing cluster
wait_for_cluster_integration() {
    local max_attempts=60  # 10 minutes
    local attempt=0
    
    echo "Waiting for cluster integration..."
    
    # Wait for ceph.conf
    while [[ ! -f /etc/ceph/ceph.conf ]]; do
        attempt=$((attempt + 1))
        if [[ $attempt -ge $max_attempts ]]; then
            echo "Error: Cluster config not deployed after $max_attempts attempts"
            return 1
        fi
        echo "Waiting for /etc/ceph/ceph.conf... (Attempt $attempt/$max_attempts)"
        sleep 10
    done
    
    # Wait for admin keyring
    attempt=0
    while [[ ! -f /etc/ceph/ceph.client.admin.keyring ]]; do
        attempt=$((attempt + 1))
        if [[ $attempt -ge $max_attempts ]]; then
            echo "Error: Admin keyring not deployed after $max_attempts attempts"
            return 1
        fi
        echo "Waiting for /etc/ceph/ceph.client.admin.keyring... (Attempt $attempt/$max_attempts)"
        sleep 10
    done
    
    # Wait for cephadm to be able to connect to cluster
    attempt=0
    while ! /usr/local/bin/cephadm shell -- ceph -s &>/dev/null; do
        attempt=$((attempt + 1))
        if [[ $attempt -ge $max_attempts ]]; then
            echo "Error: Cannot connect to cluster after $max_attempts attempts"
            return 1
        fi
        echo "Waiting for cluster connectivity... (Attempt $attempt/$max_attempts)"
        sleep 10
    done
    
    echo "Node successfully integrated into cluster"
    return 0
}

echo "Running in '$MODE' mode..."

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

# Handle different modes
if [[ "$MODE" == "first" ]]; then
    echo "First node mode: Bootstrapping new Ceph cluster..."
    
    # Check if this node already has Ceph installed
    if [[ -f /etc/ceph/ceph.conf ]]; then
        echo "Ceph already configured on this node"
        
        # Check if OSDs already exist on this node using Ceph commands
        HOSTNAME="$(hostname)"
        existing_osds_count=0
        
        # Try to get OSD count for this host from Ceph
        if existing_osds_output=$(/usr/local/bin/cephadm shell -- ceph orch ps --daemon_type osd --hostname "${HOSTNAME}" --format json 2>/dev/null); then
            existing_osds_count=$(echo "${existing_osds_output}" | python3 -c "
import json, sys
try:
    data = json.load(sys.stdin)
    print(len(data) if isinstance(data, list) else 0)
except:
    print(0)
" 2>/dev/null || echo "0")
        fi
        
        if [[ ${existing_osds_count} -gt 0 ]]; then
            echo "Found ${existing_osds_count} existing OSDs on this node. Nothing to do."
            echo "Ceph cluster status:"
            /usr/local/bin/cephadm shell -- ceph -s 2>/dev/null || echo "Could not get cluster status"
            exit 0
        else
            echo "No OSDs found on this node (checked via ceph orch ps) - will create OSDs only"
            SKIP_BOOTSTRAP=true
            SKIP_OSD_CREATION=false
        fi
    else
        echo "Bootstrapping new cluster with IP: ${INTERNAL_IP}"
        MON_IP="${INTERNAL_IP}"
        SKIP_BOOTSTRAP=false
        SKIP_OSD_CREATION=false
    fi
    
elif [[ "$MODE" == "others" ]]; then
    echo "Others mode: Joining existing Ceph cluster..."
    
    # Check if this node already has Ceph installed
    if [[ -f /etc/ceph/ceph.conf ]]; then
        echo "Ceph already configured on this node"
        
        # Check if OSDs already exist on this node using Ceph commands
        HOSTNAME="$(hostname)"
        existing_osds_count=0
        
        # Try to get OSD count for this host from Ceph
        if existing_osds_output=$(/usr/local/bin/cephadm shell -- ceph orch ps --daemon_type osd --hostname "${HOSTNAME}" --format json 2>/dev/null); then
            existing_osds_count=$(echo "${existing_osds_output}" | python3 -c "
import json, sys
try:
    data = json.load(sys.stdin)
    print(len(data) if isinstance(data, list) else 0)
except:
    print(0)
" 2>/dev/null || echo "0")
        fi
        
        if [[ ${existing_osds_count} -gt 0 ]]; then
            echo "Found ${existing_osds_count} existing OSDs on this node. Nothing to do."
            exit 0
        else
            echo "No OSDs found on this node (checked via ceph orch ps) - will create OSDs only"
            SKIP_BOOTSTRAP=true
            SKIP_OSD_CREATION=false
        fi
    else
        echo "This node needs to be added to the cluster by the orchestrator."
        echo "After running this script on the first node, add this node using:"
        echo "  ceph orch host add $(hostname) $(hostname -I | awk '{print $1}') --labels _admin"
        echo ""
        echo "Then wait for the orchestrator to deploy configuration files."
        echo "You can re-run this script after the node is added to create OSDs."
        
        # Exit here - don't try to create OSDs yet
        exit 0
    fi
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
    echo "Skipping bootstrap - node prepared for cluster joining"
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
    # For nodes that were just added, wait for integration
    if [[ "$MODE" == "others" ]] && [[ ! -f /etc/ceph/ceph.conf ]]; then
        echo "Waiting for orchestrator to deploy cluster configuration..."
        if ! wait_for_cluster_integration; then
            echo "Failed to integrate with cluster. Please check orchestrator status."
            exit 1
        fi
    fi
    
    # Check if we have cluster config (required for OSD creation)
    if [[ ! -f /etc/ceph/ceph.conf ]]; then
        echo "No cluster configuration found. Node needs to be added to cluster first."
        echo "OSDs will be created when this script runs again after cluster joining."
        exit 0
    fi
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
    echo "DEBUG: Processing device detection results:"
    while IFS=':' read -r device rotational; do
        if [[ -n "${device}" ]]; then
            echo "DEBUG: Device ${device} -> rotational=${rotational}"
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
    echo "No SSD devices found. Cannot create OSDs without SSDs for WAL/DB."
    echo "Available HDDs: ${#HDD_DEVICES[@]}"
    echo "Required: At least 1 SSD for every ${HDDS_PER_SSD} HDDs"
    echo "Skipping OSD creation to respect HDD per SSD ratio requirements."
    exit 0
else
    echo "Creating OSDs with shared DB configuration (up to ${HDDS_PER_SSD} HDDs per SSD)..."
    
    # Calculate maximum HDDs we can use based on available SSDs
    max_usable_hdds=$((${#SSD_DEVICES[@]} * HDDS_PER_SSD))
    actual_hdds_to_use=${#HDD_DEVICES[@]}
    
    if [[ ${#HDD_DEVICES[@]} -gt ${max_usable_hdds} ]]; then
        actual_hdds_to_use=${max_usable_hdds}
        echo "Found ${#HDD_DEVICES[@]} HDDs but can only use ${actual_hdds_to_use} HDDs with ${#SSD_DEVICES[@]} SSDs (${HDDS_PER_SSD} HDDs per SSD)."
        excess_hdd_count=$((${#HDD_DEVICES[@]} - actual_hdds_to_use))
        echo "Ignoring ${excess_hdd_count} excess HDDs: ${HDD_DEVICES[@]:${actual_hdds_to_use}}"
    fi
    
    # Process HDDs in groups, assigning one SSD for DB per group
    hdd_index=0
    ssd_index=0
    
    while [[ ${hdd_index} -lt ${actual_hdds_to_use} && ${ssd_index} -lt ${#SSD_DEVICES[@]} ]]; do
        # Calculate how many HDDs to assign to this SSD (respect HDDS_PER_SSD limit)
        remaining_hdds=$((actual_hdds_to_use - hdd_index))
        hdds_for_this_ssd=$((remaining_hdds < HDDS_PER_SSD ? remaining_hdds : HDDS_PER_SSD))
        
        # Build comma-separated list of HDDs for this batch
        data_devices_list=""
        for ((i=0; i<hdds_for_this_ssd; i++)); do
            if [[ $((hdd_index + i)) -lt ${actual_hdds_to_use} ]]; then
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
    
    # Report any unused SSDs
    if [[ ${ssd_index} -lt ${#SSD_DEVICES[@]} ]]; then
        unused_ssd_count=$((${#SSD_DEVICES[@]} - ssd_index))
        echo "Note: ${unused_ssd_count} SSDs remain unused: ${SSD_DEVICES[@]:${ssd_index}}"
        echo "These SSDs could be used for additional HDD groups if more HDDs were available."
    fi
fi

echo "OSD creation complete on ${HOSTNAME}. Cluster status:"
/usr/local/bin/cephadm shell -- ceph -s
