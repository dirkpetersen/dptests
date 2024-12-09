#!/bin/bash

# Function to check if service is already installed
check_service_installed() {
    local service_name="aws-sso-cache-refresh"
    if systemctl --user list-unit-files | grep -q "^$service_name.service"; then
        return 0
    fi
    return 1
}

# Function to check if running as a systemd service
check_systemd() {
    if [[ -n "${INVOCATION_ID:-}" ]] || [[ -n "${USER_INVOCATION_ID:-}" ]]; then
        return 0
    fi
    return 1
}

# Function to remove systemd user service
remove_service() {
    local service_name="aws-sso-cache-refresh"
    
    # Disable and stop the service
    systemctl --user disable "$service_name.service" 2>/dev/null
    systemctl --user stop "$service_name.service" 2>/dev/null
    
    # Remove service file
    rm -f "$HOME/.config/systemd/user/$service_name.service"
    
    # Remove from crontab
    crontab -l 2>/dev/null | grep -v "$service_name" | crontab -
    
    # Reload systemd daemon
    systemctl --user daemon-reload
    
    echo "Service removed successfully!"
}

# Function to install as systemd user service
install_service() {
    local script_path="$1"
    local service_name="aws-sso-cache-refresh"
    local service_dir="$HOME/.config/systemd/user"
    
    # Create systemd user directory if it doesn't exist
    mkdir -p "$service_dir"
    
    # Create service file
    cat > "$service_dir/$service_name.service" <<EOF
[Unit]
Description=AWS SSO Cache Refresh Service
After=network.target

[Service]
Type=simple
Environment="PATH=/usr/local/bin:/usr/bin:/bin"
Environment="HOME=$HOME"
ExecStart=/bin/bash -c "$script_path"
Restart=always
RestartSec=60
StandardOutput=journal
StandardError=journal
RemainAfterExit=no

[Install]
WantedBy=default.target
EOF

    # Reload systemd daemon
    systemctl --user daemon-reload

    # Enable and start the service
    systemctl --user enable "$service_name.service"
    systemctl --user restart "$service_name.service"

    # Add to user crontab if not already present
    if ! crontab -l 2>/dev/null | grep -q "$service_name"; then
        (crontab -l 2>/dev/null; echo "@reboot systemctl --user start $service_name.service") | crontab -
    fi

    echo "Service installed and started successfully!"
    echo "It will also run at boot time via crontab"
    
    # Show service status
    systemctl --user status "$service_name.service"
}

# Function to check if jq is installed
check_jq() {
    if ! command -v jq &> /dev/null
    then
        echo "Error: jq is not installed. Please install it first."
        echo "On Ubuntu/Debian: sudo apt-get install jq"
        echo "On MacOS: brew install jq"
        exit 1
    fi
}

# Function to find the most recent SSO cache file
find_latest_cache_file() {
    local cache_dir="$HOME/.aws/sso/cache"
    local latest_file

    if [ ! -d "$cache_dir" ]
    then
        echo "Error: AWS SSO cache directory not found"
        exit 1
    fi

    latest_file=$(find "$cache_dir" -name "*.json" -type f -exec ls -t {} + | head -n 1)
    
    if [ -z "$latest_file" ]
    then
        echo "Error: No SSO cache files found"
        exit 1
    fi

    echo "$latest_file"
}

# Function to read and parse the cache file
parse_cache_file() {
    local cache_file="$1"
    local json_content

    if [[ ! -f "$cache_file" ]]; then 
        echo "Error: Cache file not found: $cache_file"
        exit 1
    fi

    json_content=$(cat "$cache_file")

    # Extract required values using jq
    local client_id=$(echo "$json_content" | jq -r '.clientId')
    local client_secret=$(echo "$json_content" | jq -r '.clientSecret')
    local refresh_token=$(echo "$json_content" | jq -r '.refreshToken')
    local region=$(echo "$json_content" | jq -r '.region // "us-west-2"')

    # Check if required values exist
    if [[ -z "$client_id" || -z "$client_secret"  || -z "$refresh_token" ]]; then
        echo "Error: Missing required values in cache file"
        exit 1
    fi

    # Return values in an array
    echo "$client_id|$client_secret|$refresh_token|$region"
}

# Function to refresh the token
refresh_token() {
    local client_id="$1"
    local client_secret="$2"
    local refresh_token="$3"
    local region="$4"

    echo "Refreshing AWS SSO token..."
    
    aws sso-oidc create-token \
        --client-id "$client_id" \
        --client-secret "$client_secret" \
        --grant-type refresh_token \
        --refresh-token "$refresh_token" \
        --region "$region"

    if [[ $? -eq 0 ]]; then
        echo "Token refresh successful!"
    else
        echo "Error: Failed to refresh token"
        exit 1
    fi
}

# Main execution
main() {
    # If --remove flag is present, remove the service and exit
    if [[ "$1" == "--remove" ]]; then
        remove_service
        exit 0
    fi

    # If not running as service, check if already installed
    if ! check_systemd; then
        if check_service_installed; then
            echo "Service is already installed!"
            echo "Use --remove option to uninstall if needed"
            exit 1
        else
            echo "Installing as systemd user service..."
            install_service "$(realpath "$0")"
            exit 0
        fi
    fi

    # Check for jq
    check_jq

    while true; do
        # Find the latest cache file
        cache_file=$(find_latest_cache_file)
        echo "Using cache file: $cache_file"

        # Parse the cache file
        IFS='|' read -r client_id client_secret refresh_token region <<< "$(parse_cache_file "$cache_file")"

        # Refresh the token
        refresh_token "$client_id" "$client_secret" "$refresh_token" "$region"

        # Wait for 30 minutes before next refresh
        echo "Waiting 1800 seconds before next refresh..."
        sleep 1800
    done
}

# Run the script
main "$@"
