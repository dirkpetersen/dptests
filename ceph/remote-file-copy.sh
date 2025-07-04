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
# Note: EC2_KEY_FILE should be set by the calling script, no default provided
: "${SSH_OPTS:="-o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null -o ConnectTimeout=10 -o ServerAliveInterval=5 -o ServerAliveCountMax=3"}"

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
    echo "DEBUG: Using local machine as relay (SSH keys only exist locally)"
    
    # Create temporary files for metadata and content
    local temp_file=$(mktemp)
    local temp_metadata=$(mktemp)
    
    # Cleanup function
    cleanup() {
        rm -f "$temp_file" "$temp_metadata" 2>/dev/null
    }
    trap cleanup EXIT
    
    # Step 1: Get file metadata (permissions, ownership) from source via local SSH
    echo "Getting file metadata from source via local SSH..."
    if ! timeout 30 ssh -i "${EC2_KEY_FILE}" ${SSH_OPTS} "${source_user}@${source_host}" \
        "sudo stat -c '%a %U %G' '${source_path}'" > "$temp_metadata" 2>/dev/null; then
        echo "Error: Failed to get metadata from ${source_user}@${source_host}:${source_path}"
        return 1
    fi
    
    # Read metadata
    read -r file_perms file_owner file_group < "$temp_metadata"
    echo "✓ Source file permissions: ${file_perms}, owner: ${file_owner}, group: ${file_group}"
    
    # Step 2: Copy file content from source to local temp file via SSH
    echo "Downloading file from source to local temp file..."
    if ! timeout 30 ssh -i "${EC2_KEY_FILE}" ${SSH_OPTS} "${source_user}@${source_host}" \
        "sudo cat '${source_path}'" > "$temp_file" 2>/dev/null; then
        echo "Error: Failed to download file from ${source_user}@${source_host}:${source_path}"
        return 1
    fi
    
    # Verify we got the file
    if [[ ! -s "$temp_file" ]]; then
        echo "Error: Downloaded file is empty or download failed"
        return 1
    fi
    echo "✓ Downloaded file successfully ($(wc -c < "$temp_file") bytes)"
    
    # Step 3: Create destination directory if it doesn't exist via local SSH
    local dest_dir=$(dirname "$dest_path")
    echo "Ensuring destination directory exists: ${dest_dir}"
    if ! timeout 30 ssh -i "${EC2_KEY_FILE}" ${SSH_OPTS} "${dest_user}@${dest_host}" \
        "sudo mkdir -p '${dest_dir}'" 2>/dev/null; then
        echo "Error: Failed to create directory ${dest_dir} on ${dest_user}@${dest_host}"
        return 1
    fi
    
    # Step 4: Upload file content to destination via local SSH
    echo "Uploading file from local temp to destination..."
    if ! timeout 30 bash -c "cat '$temp_file' | ssh -i '${EC2_KEY_FILE}' ${SSH_OPTS} '${dest_user}@${dest_host}' 'sudo tee \"${dest_path}\" >/dev/null'"; then
        echo "Error: Failed to upload file to ${dest_user}@${dest_host}:${dest_path}"
        return 1
    fi
    
    # Step 5: Set permissions and ownership on destination via local SSH
    echo "Setting permissions (${file_perms}) and ownership (${file_owner}:${file_group}) on destination..."
    if ! timeout 30 ssh -i "${EC2_KEY_FILE}" ${SSH_OPTS} "${dest_user}@${dest_host}" \
        "sudo chmod '${file_perms}' '${dest_path}' && sudo chown '${file_owner}:${file_group}' '${dest_path}'" 2>/dev/null; then
        echo "Warning: Failed to set permissions/ownership on ${dest_user}@${dest_host}:${dest_path}"
        echo "File copied but may have incorrect permissions"
    fi
    
    echo "✓ Successfully copied ${source} to ${destination} via local relay"
    return 0
}

# Specialized function for SSH key distribution (appending to authorized_keys)
# Enhanced version supports different login user vs target user
# Usage: remote_ssh_key_copy 'user1@server1:/path/key.pub' 'user3@server2' 'root'
# Usage: remote_ssh_key_copy 'user1@server1:/path/key.pub' 'user3@server2'  # defaults to root target
remote_ssh_key_copy() {
    local source="$1" 
    local dest_user_host="$2"
    local target_user="${3:-root}"  # Target user whose authorized_keys to update (defaults to root)
    
    if [[ -z "$source" || -z "$dest_user_host" ]]; then
        echo "Error: Usage: remote_ssh_key_copy 'user@host:/path/to/key.pub' 'login_user@dest_host' [target_user]"
        echo "  source: user@host:/path/to/key.pub - where to get the SSH public key"
        echo "  dest_user_host: login_user@dest_host - user who can login and has sudo access"
        echo "  target_user: user whose authorized_keys to update (optional, defaults to 'root')"
        echo ""
        echo "Example: remote_ssh_key_copy 'user1@server1:/etc/ceph/ceph.pub' 'user3@server2' 'root'"
        return 1
    fi
    
    echo "DEBUG: Using SSH key file: ${EC2_KEY_FILE}"
    echo "DEBUG: Using local machine as relay for SSH key copy"
    
    # Parse source
    local source_user_host="${source%:*}"
    local source_path="${source#*:}"
    local source_user="${source_user_host%@*}"
    local source_host="${source_user_host#*@}"
    
    # Parse destination (login user)
    local login_user="${dest_user_host%@*}"
    local dest_host="${dest_user_host#*@}"
    
    # Validate parsing
    if [[ "$source_user_host" == "$source_path" ]]; then
        echo "Error: Invalid source format. Use 'user@host:/path/to/key.pub'"
        return 1
    fi
    
    echo "Copying SSH key from ${source_user}@${source_host}:${source_path}"
    echo "  to ${target_user}@${dest_host}:.ssh/authorized_keys"
    echo "  via login user ${login_user}@${dest_host}"
    echo "All SSH connections originate from local machine where keys exist"
    
    # Create temporary file for the key
    local temp_key=$(mktemp)
    trap "rm -f '$temp_key'" EXIT
    
    # Step 1: Download SSH public key from source to local temp file
    echo "Downloading SSH public key from source via local SSH..."
    if ! timeout 30 ssh -i "${EC2_KEY_FILE}" ${SSH_OPTS} "${source_user}@${source_host}" \
        "sudo cat '${source_path}'" > "$temp_key" 2>/dev/null; then
        echo "Error: Failed to download SSH key from ${source_user}@${source_host}:${source_path}"
        return 1
    fi
    
    # Verify we got a valid key
    if [[ ! -s "$temp_key" ]]; then
        echo "Error: Downloaded SSH key is empty or download failed"
        return 1
    fi
    
    # Basic validation that this looks like an SSH public key
    if ! grep -qE '^(ssh-rsa|ssh-ed25519|ssh-ecdsa|ecdsa-sha2-)' "$temp_key"; then
        echo "Warning: Downloaded content doesn't appear to be a valid SSH public key"
        echo "First line: $(head -1 "$temp_key")"
    fi
    
    echo "✓ SSH key downloaded successfully ($(wc -c < "$temp_key") bytes)"
    
    # Step 2: Determine target user's home directory and authorized_keys path
    local target_home_dir
    local authorized_keys_path
    
    if [[ "$target_user" == "root" ]]; then
        target_home_dir="/root"
        authorized_keys_path="/root/.ssh/authorized_keys"
    else
        # Get target user's home directory
        echo "Getting home directory for user ${target_user}..."
        target_home_dir=$(timeout 30 ssh -i "${EC2_KEY_FILE}" ${SSH_OPTS} "${login_user}@${dest_host}" \
            "sudo getent passwd '${target_user}' | cut -d: -f6" 2>/dev/null)
        
        if [[ -z "$target_home_dir" ]]; then
            echo "Error: Could not determine home directory for user ${target_user}"
            return 1
        fi
        
        authorized_keys_path="${target_home_dir}/.ssh/authorized_keys"
        echo "✓ Target user ${target_user} home directory: ${target_home_dir}"
    fi
    
    # Step 3: Check if key already exists to avoid duplicates
    echo "Checking if SSH key already exists in ${target_user}'s authorized_keys..."
    local key_content
    key_content=$(cat "$temp_key" | tr -d '\n\r' | awk '{print $1 " " $2}')  # Get key type and key data only
    
    if timeout 30 ssh -i "${EC2_KEY_FILE}" ${SSH_OPTS} "${login_user}@${dest_host}" \
        "sudo test -f '${authorized_keys_path}' && sudo grep -Fq '${key_content}' '${authorized_keys_path}'" 2>/dev/null; then
        echo "✓ SSH key already exists in ${target_user}@${dest_host}:${authorized_keys_path}, skipping duplicate"
        return 0
    fi
    
    # Step 4: Upload SSH key to destination and append to authorized_keys
    echo "Uploading SSH key to destination and adding to ${target_user}'s authorized_keys..."
    
    # Create .ssh directory with proper permissions for target user
    if ! timeout 30 ssh -i "${EC2_KEY_FILE}" ${SSH_OPTS} "${login_user}@${dest_host}" \
        "sudo mkdir -p '${target_home_dir}/.ssh' && sudo chown '${target_user}:${target_user}' '${target_home_dir}/.ssh' && sudo chmod 700 '${target_home_dir}/.ssh'" 2>/dev/null; then
        echo "Error: Failed to create/setup .ssh directory for ${target_user}"
        return 1
    fi
    
    # Append key to authorized_keys with proper ownership and permissions
    if cat "$temp_key" | timeout 30 ssh -i "${EC2_KEY_FILE}" ${SSH_OPTS} "${login_user}@${dest_host}" \
        "sudo tee -a '${authorized_keys_path}' >/dev/null && sudo chown '${target_user}:${target_user}' '${authorized_keys_path}' && sudo chmod 600 '${authorized_keys_path}'" 2>/dev/null; then
        echo "✓ Successfully added SSH key to ${target_user}@${dest_host}:${authorized_keys_path}"
        
        # Verify the key was added
        local key_count
        key_count=$(timeout 30 ssh -i "${EC2_KEY_FILE}" ${SSH_OPTS} "${login_user}@${dest_host}" \
            "sudo wc -l < '${authorized_keys_path}'" 2>/dev/null)
        echo "✓ Total SSH keys in ${target_user}'s authorized_keys: ${key_count:-unknown}"
        
        return 0
    else
        echo "Error: Failed to add SSH key to ${target_user}@${dest_host}:${authorized_keys_path}"
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
    echo "# Copy SSH key with different login user and target user"
    echo "remote_ssh_key_copy 'user1@server1:/etc/xxx/xxx.pub' 'user3@server2' 'root'"
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
    elif [[ "$1" == "ssh-key" && $# -eq 4 ]]; then
        remote_ssh_key_copy "$2" "$3" "$4"
    else
        echo "Usage:"
        echo "  $0 copy 'user@host:/source' 'user@host:/dest'"
        echo "  $0 ssh-key 'user@host:/key.pub' 'user@host' [target_user]"
        echo ""
        example_usage
    fi
fi