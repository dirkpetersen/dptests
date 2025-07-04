#!/bin/bash

CEPH_RELEASE=19.2.2  # replace this with the active release
#CEPH_RELEASE=18.2.7

# Configuration: Number of HDDs that should share one SSD for DB
: "${HDDS_PER_SSD:=2}"

# Accept mode as first argument: 'first', 'others', or 'join_cluster'
MODE="$1"

if [[ "$MODE" != "first" && "$MODE" != "others" && "$MODE" != "join_cluster" && "$MODE" != "create_osds" ]]; then
    echo "Usage: $0 {first|others|join_cluster|create_osds}"
    echo "  first        - Bootstrap the first node (MON node) of a new Ceph cluster"
    echo "  others       - Prepare an additional node to join the Ceph cluster"
    echo "  join_cluster - Add this node to cluster via orchestrator (called from first node)"
    echo "  create_osds  - Create OSDs on target hosts (run from first node)"
    exit 1
fi

# ==============================================================================
# UTILITY FUNCTIONS
# ==============================================================================

# Function to execute ceph orch commands with debug logging
ceph_orch_debug() {
    local cmd="$*"
    echo "DEBUG: Executing ceph orch command: $cmd"
    /usr/local/bin/cephadm shell -- ceph orch $cmd
}

# Function to execute ceph commands with debug logging
ceph_debug() {
    local cmd="$*"
    echo "DEBUG: Executing ceph command: $cmd"
    /usr/local/bin/cephadm shell -- ceph $cmd
}

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

# ==============================================================================
# INITIALIZATION FUNCTIONS
# ==============================================================================

# Download and configure cephadm
setup_cephadm() {
    echo "Setting up cephadm..."
    curl -o /usr/local/bin/cephadm https://download.ceph.com/rpm-${CEPH_RELEASE}/el9/noarch/cephadm
    chmod +x /usr/local/bin/cephadm
}

# Get internal IP address
get_internal_ip() {
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
    echo "Internal IP: ${INTERNAL_IP}"
}

# Wait for podman installation
wait_for_podman() {
    echo "Waiting for podman to be available (installed by cloud-init)..."
    local podman_check_attempts=30 # Approx 5 minutes (30 * 10s)
    local podman_check_count=0
    
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
}

# Configure udev rules for AWS EBS volumes
configure_udev_rules() {
    echo "Configuring udev rules for AWS EBS volumes..."
    cat > /etc/udev/rules.d/99-aws-ebs.rules << 'EOF'
# Add this content (targets ONLY "Amazon Elastic Block Store"):
SUBSYSTEM=="block", KERNEL=="nvme*", ENV{ID_MODEL}=="Amazon Elastic Block Store", ENV{DEVTYPE}=="disk", ATTR{queue/rotational}="1"
EOF

    # Reload udev rules
    echo "Reloading udev rules..."
    udevadm control --reload-rules
    udevadm trigger
}

# Source remote file copy functions
source_remote_functions() {
    if [[ -f "$(dirname "$0")/remote-file-copy.sh" ]]; then
        source "$(dirname "$0")/remote-file-copy.sh"
    elif [[ -f "./remote-file-copy.sh" ]]; then
        source "./remote-file-copy.sh"
    fi
}

# ==============================================================================
# MODE-SPECIFIC FUNCTIONS
# ==============================================================================

# Handle 'first' mode - bootstrap new cluster
handle_first_mode() {
    echo "First node mode: Bootstrapping new Ceph cluster..."
    
    # Check if this node already has Ceph installed
    if [[ -f /etc/ceph/ceph.conf ]]; then
        echo "Ceph already configured on this node"
        
        # Check if OSDs already exist on this node using Ceph commands
        HOSTNAME="$(hostname)"
        existing_osds_count=0
        
        # Try to get OSD count for this host from Ceph
        if existing_osds_output=$(ceph_orch_debug "ps --daemon_type osd --hostname ${HOSTNAME} --format json" 2>/dev/null); then
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
}

# Handle 'others' mode - prepare node to join cluster
handle_others_mode() {
    echo "Others mode: Preparing node to join Ceph cluster..."
    
    # This mode just prepares the node - it doesn't create OSDs
    # OSD creation is handled by the orchestrator from the first node
    
    HOSTNAME="$(hostname)"
    echo "########### Node ${HOSTNAME} prepared and ready to join cluster."
    echo "This node should be added to cluster from the first node using:"
    echo "  ceph orch host add ${HOSTNAME} ${INTERNAL_IP}"
    echo ""
    echo "After adding to cluster, OSDs will be created from the first node using:"
    echo "  ceph orch daemon add osd ${HOSTNAME}:data_devices=<devices>,db_devices=<ssd>,osds_per_device=1"
    
    # Always skip bootstrap and OSD creation in 'others' mode
    SKIP_BOOTSTRAP=true
    SKIP_OSD_CREATION=true
}

# Handle 'join_cluster' mode - add node to existing cluster
handle_join_cluster_mode() {
    echo "Join cluster mode: Adding node to existing cluster via orchestrator..."
    
    # Get parameters from environment variables (set by calling script)
    : "${FIRST_NODE_INTERNAL_IP:?}"
    : "${TARGET_HOSTNAME:?}"
    : "${TARGET_INTERNAL_IP:?}"
    
    echo "Adding host ${TARGET_HOSTNAME} (${TARGET_INTERNAL_IP}) to cluster..."
    echo "First node: ${FIRST_NODE_INTERNAL_IP}"
    
    # SSH key distribution is handled by the local machine before this runs
    echo "SSH key distribution to ${TARGET_HOSTNAME} is handled by local machine prior to join_cluster mode."
    echo "Proceeding with host addition to cluster..."
    
    # Add the host to the cluster using ceph orchestrator
    echo "DEBUG: Executing ceph orch command: host add ${TARGET_HOSTNAME} ${TARGET_INTERNAL_IP}"
    if /usr/local/bin/cephadm shell -- ceph orch host add ${TARGET_HOSTNAME} ${TARGET_INTERNAL_IP}; then
        echo "Successfully added ${TARGET_HOSTNAME} to cluster."
        
        # Wait for the node to be fully joined to the cluster
        wait_for_host_integration "${TARGET_HOSTNAME}"
        
        echo "Node successfully joined cluster."
        exit 0
    else
        echo "Error: Failed to add ${TARGET_HOSTNAME} to cluster."
        echo "Manual addition required: ceph orch host add ${TARGET_HOSTNAME} ${TARGET_INTERNAL_IP}"
        exit 1
    fi
}

# Handle 'create_osds' mode - create OSDs on target host
handle_create_osds_mode() {
    echo "Create OSDs mode: Creating OSDs on target host from orchestrator..."
    
    # Get parameters from environment variables (set by calling script)
    : "${TARGET_HOSTNAME:?}"
    
    echo "########### Creating OSDs on ${TARGET_HOSTNAME} from orchestrator..."
    
    # This mode runs on the first node and creates OSDs on the target host
    SKIP_BOOTSTRAP=true
    SKIP_OSD_CREATION=false
}

# Wait for host to be integrated into cluster
wait_for_host_integration() {
    local hostname="$1"
    echo "Waiting for ${hostname} to be fully integrated into the cluster..."
    
    local join_attempts=30
    local join_count=0
    
    while [[ ${join_count} -lt ${join_attempts} ]]; do
        # Check if the node appears in host list
        echo "DEBUG: Executing ceph orch command: host ls --format json"
        if host_list_output=$(/usr/local/bin/cephadm shell -- ceph orch host ls --format json 2>/dev/null); then
            if echo "${host_list_output}" | jq -r '.[].hostname' 2>/dev/null | grep -q "^${hostname}$"; then
                echo "${hostname} is now fully integrated into the cluster."
                break
            fi
        fi
        
        join_count=$((join_count + 1))
        echo "Waiting for ${hostname} to appear in cluster host list... (Attempt ${join_count}/${join_attempts})"
        sleep 10
    done
    
    if [[ ${join_count} -eq ${join_attempts} ]]; then
        echo "Warning: ${hostname} did not fully join the cluster within expected time."
        exit 1
    fi
}

# ==============================================================================
# BOOTSTRAP FUNCTIONS
# ==============================================================================

# Bootstrap new Ceph cluster
bootstrap_cluster() {
    if [[ "${SKIP_BOOTSTRAP}" == "false" ]]; then
        echo "Bootstrapping Ceph cluster with MON IP: ${MON_IP}"
        /usr/local/bin/cephadm bootstrap --mon-ip "${MON_IP}" --ssh-user root

        # Restart ceph services after bootstrap
        echo "Restarting Ceph services..."
        systemctl restart ceph.target

        # Wait for Ceph cluster to be ready
        echo "Waiting for Ceph cluster to be ready..."
        sleep 30
        
        # Ensure SSH keys are generated for orchestrator
        ensure_ssh_keys
    else
        echo "Skipping bootstrap - node prepared for cluster joining"
    fi
}

# Ensure SSH keys are available for orchestrator
ensure_ssh_keys() {
    echo "Ensuring SSH keys are available for orchestrator..."
    if [[ ! -f /etc/ceph/ceph.pub ]]; then
        echo "Generating SSH keys for Ceph orchestrator..."
        echo "DEBUG: Executing ceph command: cephadm generate-key"
        /usr/local/bin/cephadm shell -- ceph cephadm generate-key || true
        
        # Alternative method: use ceph orch to generate keys
        echo "DEBUG: Executing ceph orch command: set backend cephadm"
        /usr/local/bin/cephadm shell -- ceph orch set backend cephadm || true
        
        # Check if key was created
        if [[ -f /etc/ceph/ceph.pub ]]; then
            echo "✓ SSH key generated successfully"
        else
            echo "⚠ SSH key generation may have failed. Manual key setup might be required."
        fi
    else
        echo "✓ SSH key already exists"
    fi
}

# ==============================================================================
# OSD CREATION FUNCTIONS
# ==============================================================================

# Check if cephadm can connect to cluster
check_cluster_connection() {
    echo "Checking if cephadm can connect to cluster..."
    
    # Try to connect to cluster using cephadm shell
    if /usr/local/bin/cephadm shell -- ceph -s >/dev/null 2>&1; then
        echo "✓ Successfully connected to Ceph cluster via cephadm"
        return 0
    else
        echo "Cannot connect to Ceph cluster. Checking possible reasons..."
        
        # Check if this host is known to the cluster
        echo "Checking if this host is added to the cluster..."
        HOSTNAME="$(hostname)"
        
        # This might fail if not added to cluster yet - that's okay
        if /usr/local/bin/cephadm shell -- ceph orch host ls --format json 2>/dev/null | jq -r '.[].hostname' 2>/dev/null | grep -q "^${HOSTNAME}$"; then
            echo "✓ Host ${HOSTNAME} is known to the cluster"
            echo "Waiting for cephadm connection to become available..."
            
            # Wait a bit for cephadm to establish connection
            sleep 30
            
            # Try connection again
            if /usr/local/bin/cephadm shell -- ceph -s >/dev/null 2>&1; then
                echo "✓ Connection established after wait"
                return 0
            else
                echo "Warning: Still cannot connect to cluster. Proceeding anyway - cephadm may establish connection during OSD creation."
                return 0
            fi
        else
            echo "Host ${HOSTNAME} is not added to the cluster yet."
            echo "This host needs to be added via: ceph orch host add ${HOSTNAME} ${INTERNAL_IP}"
            echo "Skipping OSD creation until host is added to cluster."
            exit 0
        fi
    fi
}

# Detect available storage devices with retry mechanism
detect_storage_devices() {
    echo "Detecting available storage devices using ceph orch device ls..."
    SSD_DEVICES=()
    HDD_DEVICES=()

    # Get cluster info
    FSID="$(grep fsid /etc/ceph/ceph.conf | awk '{print $3}')"
    
    # Use target hostname if in create_osds mode, otherwise use current hostname
    if [[ "$MODE" == "create_osds" ]]; then
        HOSTNAME="${TARGET_HOSTNAME}"
        echo "Creating OSDs for target host: ${HOSTNAME}"
    else
        HOSTNAME="$(hostname)"
    fi

    # Retry mechanism for device detection
    local max_device_attempts=12  # 12 attempts * 10s = 2 minutes
    local device_attempt=0
    local devices_found=false
    
    while [[ ${device_attempt} -lt ${max_device_attempts} && "${devices_found}" == "false" ]]; do
        device_attempt=$((device_attempt + 1))
        
        echo "DEBUG: Device detection attempt ${device_attempt}/${max_device_attempts}"
        echo "DEBUG: Executing ceph orch command: device ls --format json"
        
        local device_output
        device_output=$(FSID="${FSID}" /usr/local/bin/cephadm shell --fsid "${FSID}" -c /etc/ceph/ceph.conf -k /etc/ceph/ceph.client.admin.keyring -- \
            ceph orch device ls ${HOSTNAME} --format json 2>/dev/null)

        if [[ $? -eq 0 && -n "${device_output}" ]]; then
            echo "DEBUG: Raw device output (attempt ${device_attempt}):"
            echo "${device_output}" | head -5

            # Analyze device status in the response
            local device_analysis
            device_analysis=$(echo "${device_output}" | python3 -c "
import json, sys
try:
    data = json.load(sys.stdin)
    total_devices = 0
    available_devices = 0
    in_use_devices = 0
    rejected_devices = 0
    device_details = []
    
    if isinstance(data, list):
        for host_data in data:
            if isinstance(host_data, dict) and 'devices' in host_data:
                for device in host_data['devices']:
                    total_devices += 1
                    path = device.get('path', 'unknown')
                    available = device.get('available', False)
                    rejected_reasons = device.get('rejected_reasons', [])
                    
                    if available and not rejected_reasons:
                        available_devices += 1
                        device_details.append(f'{path}:available')
                    elif rejected_reasons:
                        rejected_devices += 1
                        # Get first rejection reason for summary
                        reason = rejected_reasons[0] if rejected_reasons else 'unknown'
                        device_details.append(f'{path}:rejected:{reason}')
                    else:
                        in_use_devices += 1
                        device_details.append(f'{path}:in_use')
            elif isinstance(host_data, dict) and 'path' in host_data:
                total_devices += 1
                path = host_data.get('path', 'unknown')
                available = host_data.get('available', False)
                rejected_reasons = host_data.get('rejected_reasons', [])
                
                if available and not rejected_reasons:
                    available_devices += 1
                    device_details.append(f'{path}:available')
                elif rejected_reasons:
                    rejected_devices += 1
                    reason = rejected_reasons[0] if rejected_reasons else 'unknown'
                    device_details.append(f'{path}:rejected:{reason}')
                else:
                    in_use_devices += 1
                    device_details.append(f'{path}:in_use')
    
    print(f'{total_devices}:{available_devices}:{in_use_devices}:{rejected_devices}')
    for detail in device_details:
        print(detail, file=sys.stderr)
        
except Exception as e:
    print('0:0:0:0')
    print(f'Error: {e}', file=sys.stderr)
" 2>/tmp/device_details.txt || echo "0:0:0:0")

            # Parse device analysis
            IFS=':' read -r device_count available_count in_use_count rejected_count <<< "${device_analysis}"
            
            echo "DEBUG: Device analysis (attempt ${device_attempt}):"
            echo "  Total devices: ${device_count}"
            echo "  Available: ${available_count}"  
            echo "  In use: ${in_use_count}"
            echo "  Rejected: ${rejected_count}"
            
            # Show device details if available
            if [[ -f /tmp/device_details.txt ]]; then
                echo "DEBUG: Device details:"
                while IFS=':' read -r dev_path dev_status dev_reason; do
                    if [[ -n "${dev_path}" ]]; then
                        case "${dev_status}" in
                            "available")
                                echo "  ${dev_path}: Available for OSD creation"
                                ;;
                            "in_use")
                                echo "  ${dev_path}: Already in use (likely by existing OSD)"
                                ;;
                            "rejected")
                                echo "  ${dev_path}: Rejected (${dev_reason})"
                                ;;
                            *)
                                echo "  ${dev_path}: Status unknown (${dev_status})"
                                ;;
                        esac
                    fi
                done < /tmp/device_details.txt
                rm -f /tmp/device_details.txt
            fi
            
            # Decide what to do based on device analysis
            if [[ ${available_count} -gt 0 ]]; then
                echo "Found ${available_count} available devices, proceeding to process them..."
                
                local available_devices
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
                echo "DEBUG: Processing available devices for OSD creation:"
                local temp_ssd_count=0
                local temp_hdd_count=0
                
                while IFS=':' read -r device rotational; do
                    if [[ -n "${device}" ]]; then
                        echo "DEBUG: Device ${device} -> rotational=${rotational}"
                        if [[ "${rotational}" == "False" ]]; then
                            SSD_DEVICES+=("${device}")
                            temp_ssd_count=$((temp_ssd_count + 1))
                            echo "Found available SSD: ${device}"
                        else
                            HDD_DEVICES+=("${device}") 
                            temp_hdd_count=$((temp_hdd_count + 1))
                            echo "Found available HDD: ${device}"
                        fi
                    fi
                done <<< "${available_devices}"
                
                echo "✓ Device detection successful: ${temp_ssd_count} SSDs, ${temp_hdd_count} HDDs"
                devices_found=true
                break
                
            elif [[ ${device_count} -gt 0 ]]; then
                # We have devices but none are available - analyze why
                if [[ ${in_use_count} -gt 0 && ${device_attempt} -ge 3 ]]; then
                    # If most devices are in use and we've tried a few times, this might be normal
                    echo "Found ${device_count} devices on ${HOSTNAME}:"
                    echo "  - ${available_count} available"
                    echo "  - ${in_use_count} in use (likely by existing OSDs)"
                    echo "  - ${rejected_count} rejected"
                    echo ""
                    echo "All devices appear to already be configured for Ceph OSDs."
                    echo "This is normal if OSDs were previously created on this host."
                    echo "No new OSDs need to be created."
                    devices_found=true  # Consider this a success - no work needed
                    break
                elif [[ ${rejected_count} -gt 0 && ${available_count} -eq 0 && ${in_use_count} -eq 0 ]]; then
                    # All devices are rejected - no point in retrying as this won't change
                    echo "Found ${device_count} devices on ${HOSTNAME} but all are rejected:"
                    echo "  - Available: ${available_count}"
                    echo "  - In use: ${in_use_count}"
                    echo "  - Rejected: ${rejected_count}"
                    echo ""
                    echo "All devices are rejected and cannot be used for OSDs."
                    echo "Common rejection reasons: existing partitions, too small, hardware issues."
                    echo "No retries needed - device conditions are permanent."
                    devices_found=true  # Consider this final - don't retry rejected devices
                    break
                elif [[ ${rejected_count} -gt 0 && $((rejected_count + in_use_count)) -eq ${device_count} && ${device_attempt} -ge 2 ]]; then
                    # Mix of rejected and in-use devices, no available devices - likely permanent condition
                    echo "Found ${device_count} devices on ${HOSTNAME}:"
                    echo "  - Available: ${available_count}"
                    echo "  - In use: ${in_use_count} (likely by existing OSDs)"
                    echo "  - Rejected: ${rejected_count} (permanent issues)"
                    echo ""
                    echo "All devices are either in use or rejected - no devices available for new OSDs."
                    echo "This is a permanent condition, no retries needed."
                    devices_found=true  # Consider this final
                    break
                elif [[ ${device_attempt} -lt 6 ]]; then
                    # Continue retrying for the first few attempts in case it's a timing issue
                    # But only if we have devices that might become available (not all rejected)
                    if [[ ${rejected_count} -lt ${device_count} ]]; then
                        echo "Devices exist but none available (attempt ${device_attempt}). May be timing issue, retrying..."
                    else
                        echo "All ${device_count} devices are rejected. Stopping retries."
                        devices_found=true  # Stop retrying rejected devices
                        break
                    fi
                else
                    # After several attempts, provide detailed status
                    echo "After ${device_attempt} attempts, found ${device_count} devices but none available:"
                    echo "  - Available: ${available_count}"
                    echo "  - In use: ${in_use_count}" 
                    echo "  - Rejected: ${rejected_count}"
                    echo "This may indicate all devices are already configured or have permanent issues."
                fi
            else
                # No devices found at all
                if [[ ${device_attempt} -lt 6 ]]; then
                    echo "No devices found in orchestrator inventory for ${HOSTNAME} (attempt ${device_attempt})"
                    echo "Orchestrator may still be discovering devices on newly added host"
                else
                    echo "No devices found after ${device_attempt} attempts. Host may not have storage devices."
                fi
            fi
        else
            echo "Failed to get device list from ceph orchestrator (attempt ${device_attempt})"
        fi
        
        # Wait before retry unless this was the last attempt
        if [[ ${device_attempt} -lt ${max_device_attempts} && "${devices_found}" == "false" ]]; then
            echo "Waiting 10 seconds before next device detection attempt..."
            sleep 10
        fi
    done
    
    # Final check and intelligent error handling
    if [[ "${devices_found}" == "false" ]]; then
        echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
        echo "DEVICE DETECTION SUMMARY for ${HOSTNAME}"
        echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
        echo "After ${max_device_attempts} attempts over $((max_device_attempts * 10 / 60)) minutes:"
        echo ""
        
        # Get final device status for summary
        local final_device_output
        final_device_output=$(FSID="${FSID}" /usr/local/bin/cephadm shell --fsid "${FSID}" -c /etc/ceph/ceph.conf -k /etc/ceph/ceph.client.admin.keyring -- \
            ceph orch device ls ${HOSTNAME} --format json 2>/dev/null)
        
        if [[ $? -eq 0 && -n "${final_device_output}" ]]; then
            local final_analysis
            final_analysis=$(echo "${final_device_output}" | python3 -c "
import json, sys
try:
    data = json.load(sys.stdin)
    total_devices = 0
    available_devices = 0
    in_use_devices = 0
    rejected_devices = 0
    
    if isinstance(data, list):
        for host_data in data:
            if isinstance(host_data, dict) and 'devices' in host_data:
                for device in host_data['devices']:
                    total_devices += 1
                    available = device.get('available', False)
                    rejected_reasons = device.get('rejected_reasons', [])
                    
                    if available and not rejected_reasons:
                        available_devices += 1
                    elif rejected_reasons:
                        rejected_devices += 1
                    else:
                        in_use_devices += 1
    
    print(f'{total_devices}:{available_devices}:{in_use_devices}:{rejected_devices}')
except:
    print('0:0:0:0')
" 2>/dev/null || echo "0:0:0:0")
            
            IFS=':' read -r final_total final_available final_in_use final_rejected <<< "${final_analysis}"
            
            if [[ ${final_total} -gt 0 ]]; then
                echo "Final device status:"
                echo "  • Total devices detected: ${final_total}"
                echo "  • Available for OSDs: ${final_available}"
                echo "  • Already in use: ${final_in_use}"
                echo "  • Rejected/unusable: ${final_rejected}"
                echo ""
                
                if [[ ${final_in_use} -gt 0 && ${final_available} -eq 0 ]]; then
                    echo "✓ CONCLUSION: All devices are already configured with OSDs"
                    echo "  This is normal for hosts that already have Ceph OSDs deployed."
                    echo "  No action needed - OSDs are already present and functional."
                    echo ""
                    echo "To verify existing OSDs on this host:"
                    echo "  ceph orch ps --daemon-type osd --hostname ${HOSTNAME}"
                    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
                    # Exit successfully since this is a normal condition
                    exit 0
                elif [[ ${final_rejected} -gt 0 && ${final_available} -eq 0 ]]; then
                    if [[ ${final_in_use} -gt 0 ]]; then
                        echo "✓ CONCLUSION: Mix of configured OSDs and rejected devices"
                        echo "  - ${final_in_use} devices already have OSDs (normal)"
                        echo "  - ${final_rejected} devices rejected (permanent issues)"
                        echo "  This is normal - working devices are already configured."
                        echo ""
                        echo "To verify existing OSDs on this host:"
                        echo "  ceph orch ps --daemon-type osd --hostname ${HOSTNAME}"
                        echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
                        exit 0
                    else
                        echo "✓ CONCLUSION: All devices are rejected and cannot be used"
                        echo "  Common reasons for rejection:"
                        echo "  - Devices have existing partition tables or filesystems"
                        echo "  - Devices are too small for Ceph requirements"
                        echo "  - Hardware issues or unsupported device types"
                        echo ""
                        echo "This is a permanent condition - no OSDs can be created."
                        echo "Manual intervention may be required if OSDs are expected."
                        echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
                        exit 0
                    fi
                else
                    echo "⚠ ISSUE: Devices exist but are not available for OSD creation"
                    echo "  This could indicate:"
                    echo "  - Timing issues with orchestrator device inventory"
                    echo "  - Network connectivity problems"
                    echo "  - Orchestrator cannot access or manage the devices"
                fi
            else
                echo "No devices found on ${HOSTNAME}"
                echo "This could indicate:"
                echo "  - Host has no additional storage devices"
                echo "  - Orchestrator cannot inventory devices on this host"
                echo "  - Host was recently added and needs more time"
            fi
        else
            echo "Cannot communicate with Ceph orchestrator to get device status"
            echo "This could indicate:"
            echo "  - Host ${HOSTNAME} is not properly joined to cluster"
            echo "  - Network connectivity issues"
            echo "  - Ceph cluster health problems"
        fi
        
        echo ""
        echo "Manual troubleshooting commands:"
        echo "  ceph orch device ls ${HOSTNAME}                    # Check device status"
        echo "  ceph orch host ls                                   # Verify host is in cluster"
        echo "  ceph orch ps --hostname ${HOSTNAME}                 # Check existing daemons"
        echo "  ceph -s                                             # Check cluster health"
        echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
        exit 1
    fi
}

# Validate devices and create OSDs
create_osds() {
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
        create_osds_with_shared_db
    fi
}

# Create OSDs with shared DB configuration
create_osds_with_shared_db() {
    # Calculate maximum HDDs we can use based on available SSDs
    local max_usable_hdds=$((${#SSD_DEVICES[@]} * HDDS_PER_SSD))
    local actual_hdds_to_use=${#HDD_DEVICES[@]}
    
    if [[ ${#HDD_DEVICES[@]} -gt ${max_usable_hdds} ]]; then
        actual_hdds_to_use=${max_usable_hdds}
        echo "Found ${#HDD_DEVICES[@]} HDDs but can only use ${actual_hdds_to_use} HDDs with ${#SSD_DEVICES[@]} SSDs (${HDDS_PER_SSD} HDDs per SSD)."
        local excess_hdd_count=$((${#HDD_DEVICES[@]} - actual_hdds_to_use))
        echo "Ignoring ${excess_hdd_count} excess HDDs: ${HDD_DEVICES[@]:${actual_hdds_to_use}}"
    fi
    
    # Process HDDs in groups, assigning one SSD for DB per group
    local hdd_index=0
    local ssd_index=0
    
    while [[ ${hdd_index} -lt ${actual_hdds_to_use} && ${ssd_index} -lt ${#SSD_DEVICES[@]} ]]; do
        # Calculate how many HDDs to assign to this SSD (respect HDDS_PER_SSD limit)
        local remaining_hdds=$((actual_hdds_to_use - hdd_index))
        local hdds_for_this_ssd=$((remaining_hdds < HDDS_PER_SSD ? remaining_hdds : HDDS_PER_SSD))
        
        # Build comma-separated list of HDDs for this batch
        local data_devices_list=""
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
            local ssd_device="${SSD_DEVICES[$ssd_index]}"
            echo "Creating $hdds_for_this_ssd OSDs with data on [$data_devices_list], DB on $ssd_device..."
            
            echo "DEBUG: Executing ceph orch command: daemon add osd $HOSTNAME:data_devices=$data_devices_list,db_devices=$ssd_device,osds_per_device=1"
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
        local unused_ssd_count=$((${#SSD_DEVICES[@]} - ssd_index))
        echo "Note: ${unused_ssd_count} SSDs remain unused: ${SSD_DEVICES[@]:${ssd_index}}"
        echo "These SSDs could be used for additional HDD groups if more HDDs were available."
    fi
}

# Show final cluster status
show_cluster_status() {
    echo "OSD creation complete on ${HOSTNAME}. Cluster status:"
    echo "DEBUG: Executing ceph command: -s"
    /usr/local/bin/cephadm shell -- ceph -s
}

# ==============================================================================
# MAIN EXECUTION FLOW
# ==============================================================================

echo "########## Running in '$MODE' mode..."

# Initialize common components
source_remote_functions
setup_cephadm
get_internal_ip
wait_for_podman
configure_udev_rules

# Handle different modes
case "$MODE" in
    "first")
        handle_first_mode
        ;;
    "others")
        handle_others_mode
        ;;
    "join_cluster")
        handle_join_cluster_mode
        ;;
    "create_osds")
        handle_create_osds_mode
        ;;
esac

# Bootstrap cluster if needed
bootstrap_cluster

# Create OSDs if not skipping
if [[ "${SKIP_OSD_CREATION}" != "true" ]]; then
    echo "Preparing for OSD creation..."
    check_cluster_connection
    detect_storage_devices
    create_osds
    show_cluster_status
else
    echo "Skipping OSD creation for this mode."
    exit 0
fi