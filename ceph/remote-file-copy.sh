#!/bin/bash

# Remote File Copy Function
# 
# This function copies files between remote hosts through the local machine,
# preserving permissions and creating necessary directories.
# 
# Usage: remote_file_copy "source_user@source_host:/path/to/source" "dest_user@dest_host:/path/to/dest"
# 
# Features:
# - Handles files requiring sudo access
# - Preserves file permissions and ownership
# - Creates target directories if they don't exist
# - Works with files under /etc, /root, and other protected directories
# - Uses SSH key authentication with configurable options

# Configuration variables (can be overridden by environment)
: "${EC2_KEY_FILE:="~/.ssh/id_rsa"}"
: "${SSH_OPTS:="-o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null"}"

remote_file_copy() {
    local source="$1"
    local destination="$2"
    
    if [[ -z "$source" || -z "$destination" ]]; then
        echo "Error: Usage: remote_file_copy 'user@host:/path/to/source' 'user@host:/path/to/dest'"
        return 1
    fi
    
    # Parse source
    local source_user_host="${source%:*}"
    local source_path="${source#*:}"
    local source_user="${source_user_host%@*}"
    local source_host="${source_user_host#*@}"
    
    # Parse destination  
    local dest_user_host="${destination%:*}"
    local dest_path="${destination#*:}"
    local dest_user="${dest_user_host%@*}"
    local dest_host="${dest_user_host#*@}"
    
    # Validate parsing
    if [[ "$source_user_host" == "$source_path" ]] || [[ "$dest_user_host" == "$dest_path" ]]; then
        echo "Error: Invalid format. Use 'user@host:/path/to/file'"
        return 1
    fi
    
    echo "Copying file from ${source_user}@${source_host}:${source_path} to ${dest_user}@${dest_host}:${dest_path}"
    
    # Create temporary files for metadata and content
    local temp_file=$(mktemp)
    local temp_metadata=$(mktemp)
    
    # Cleanup function
    cleanup() {
        rm -f "$temp_file" "$temp_metadata" 2>/dev/null
    }
    trap cleanup EXIT
    
    # Step 1: Get file metadata (permissions, ownership) from source
    echo "Getting file metadata from source..."
    if ! ssh -i "${EC2_KEY_FILE}" ${SSH_OPTS} "${source_user}@${source_host}" \
        "sudo stat -c '%a %U %G' '${source_path}'" > "$temp_metadata" 2>/dev/null; then
        echo "Error: Failed to get metadata from ${source_user}@${source_host}:${source_path}"
        return 1
    fi
    
    # Read metadata
    read -r file_perms file_owner file_group < "$temp_metadata"
    echo "Source file permissions: ${file_perms}, owner: ${file_owner}, group: ${file_group}"
    
    # Step 2: Copy file content from source to local temp file
    echo "Copying file content from source..."
    if ! ssh -i "${EC2_KEY_FILE}" ${SSH_OPTS} "${source_user}@${source_host}" \
        "sudo cat '${source_path}'" > "$temp_file" 2>/dev/null; then
        echo "Error: Failed to read file from ${source_user}@${source_host}:${source_path}"
        return 1
    fi
    
    # Step 3: Create destination directory if it doesn't exist
    local dest_dir=$(dirname "$dest_path")
    echo "Ensuring destination directory exists: ${dest_dir}"
    if ! ssh -i "${EC2_KEY_FILE}" ${SSH_OPTS} "${dest_user}@${dest_host}" \
        "sudo mkdir -p '${dest_dir}'" 2>/dev/null; then
        echo "Error: Failed to create directory ${dest_dir} on ${dest_user}@${dest_host}"
        return 1
    fi
    
    # Step 4: Copy file content to destination
    echo "Copying file content to destination..."
    if ! cat "$temp_file" | ssh -i "${EC2_KEY_FILE}" ${SSH_OPTS} "${dest_user}@${dest_host}" \
        "sudo tee '${dest_path}' >/dev/null"; then
        echo "Error: Failed to write file to ${dest_user}@${dest_host}:${dest_path}"
        return 1
    fi
    
    # Step 5: Set permissions and ownership on destination
    echo "Setting permissions (${file_perms}) and ownership (${file_owner}:${file_group}) on destination..."
    if ! ssh -i "${EC2_KEY_FILE}" ${SSH_OPTS} "${dest_user}@${dest_host}" \
        "sudo chmod '${file_perms}' '${dest_path}' && sudo chown '${file_owner}:${file_group}' '${dest_path}'" 2>/dev/null; then
        echo "Warning: Failed to set permissions/ownership on ${dest_user}@${dest_host}:${dest_path}"
        echo "File copied but may have incorrect permissions"
    fi
    
    echo "Successfully copied ${source} to ${destination}"
    return 0
}

# Specialized function for SSH key distribution (appending to authorized_keys)
remote_ssh_key_copy() {
    local source="$1" 
    local dest_user_host="$2"
    
    if [[ -z "$source" || -z "$dest_user_host" ]]; then
        echo "Error: Usage: remote_ssh_key_copy 'user@host:/path/to/key.pub' 'user@host'"
        return 1
    fi
    
    # Parse source
    local source_user_host="${source%:*}"
    local source_path="${source#*:}"
    local source_user="${source_user_host%@*}"
    local source_host="${source_user_host#*@}"
    
    # Parse destination
    local dest_user="${dest_user_host%@*}"
    local dest_host="${dest_user_host#*@}"
    
    echo "Copying SSH key from ${source_user}@${source_host}:${source_path} to ${dest_user}@${dest_host}:/root/.ssh/authorized_keys"
    
    # Create temporary file for the key
    local temp_key=$(mktemp)
    trap "rm -f '$temp_key'" EXIT
    
    # Step 1: Get the SSH public key from source
    echo "Reading SSH public key from source..."
    if ! ssh -i "${EC2_KEY_FILE}" ${SSH_OPTS} "${source_user}@${source_host}" \
        "sudo cat '${source_path}'" > "$temp_key" 2>/dev/null; then
        echo "Error: Failed to read SSH key from ${source_user}@${source_host}:${source_path}"
        return 1
    fi
    
    # Step 2: Set up SSH directory and append key to authorized_keys on destination
    echo "Setting up SSH directory and appending key to authorized_keys..."
    if cat "$temp_key" | ssh -i "${EC2_KEY_FILE}" ${SSH_OPTS} "${dest_user}@${dest_host}" \
        "sudo mkdir -p /root/.ssh && sudo chmod 700 /root/.ssh && sudo tee -a /root/.ssh/authorized_keys >/dev/null && sudo chmod 600 /root/.ssh/authorized_keys"; then
        echo "Successfully added SSH key to ${dest_user}@${dest_host}:/root/.ssh/authorized_keys"
        return 0
    else
        echo "Error: Failed to add SSH key to ${dest_user}@${dest_host}"
        return 1
    fi
}

# Example usage functions
example_usage() {
    echo "Example usage:"
    echo ""
    echo "# Copy a regular file"
    echo "remote_file_copy 'rocky@host1:/etc/ceph/ceph.conf' 'rocky@host2:/etc/ceph/ceph.conf'"
    echo ""
    echo "# Copy SSH key (appending to authorized_keys)"  
    echo "remote_ssh_key_copy 'rocky@host1:/etc/ceph/ceph.pub' 'rocky@host2'"
    echo ""
    echo "# Set custom SSH key file"
    echo "EC2_KEY_FILE='~/.ssh/my-key.pem' remote_file_copy 'user@host1:/path/file' 'user@host2:/path/file'"
}

# If script is executed directly, show usage
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    if [[ $# -eq 0 ]]; then
        example_usage
    elif [[ "$1" == "copy" && $# -eq 3 ]]; then
        remote_file_copy "$2" "$3"
    elif [[ "$1" == "ssh-key" && $# -eq 3 ]]; then
        remote_ssh_key_copy "$2" "$3"
    else
        echo "Usage:"
        echo "  $0 copy 'user@host:/source' 'user@host:/dest'"
        echo "  $0 ssh-key 'user@host:/key.pub' 'user@host'"
        echo ""
        example_usage
    fi
fi