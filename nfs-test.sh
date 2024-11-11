#! /bin/bash

# Script to setup NFS server or client for testing
# usage: ./nfs-test.sh [nfs_server_fqdn]
# when invoked without arguments, it sets up NFS server
# when invoked with an argument, it sets up NFS client 
# mounts server and runs scratch-dna benchmark

set -e

# Function to install common dependencies
install_common_deps() {
    sudo yum update -y
    sudo yum install -y nfs-utils rpcbind
}

# Function to install and configure NFS server
setup_nfs_server() {
    echo "Setting up NFS server..."
    
    # Install dependencies
    install_common_deps
    
    # Start and enable services
    sudo systemctl start rpcbind nfs-server
    sudo systemctl enable rpcbind nfs-server
    
    # Configure exports if directories exist
    if [ -d "/opt" ] || [ -d "/mnt/scratch" ]; then
        echo "Configuring exports..."
        
        # Create exports entries
        if [ -d "/opt" ]; then
            echo "/opt *(rw,sync,no_root_squash,no_subtree_check)" | sudo tee -a /etc/exports
            mkdir -p /opt/temp
            chown 1000:1000 /opt/temp
        fi
        
        if [ -d "/mnt/scratch" ]; then
            echo "/mnt/scratch *(rw,sync,no_root_squash,no_subtree_check)" | sudo tee -a /etc/exports
            mkdir -p /mnt/scratch/temp
            chown 1000:1000 /mnt/scratch/temp
        fi
        
        # Apply exports
        sudo exportfs -ra
    else
        echo "Neither /opt nor /mnt/scratch directories exist. Nothing to export."
        exit 1
    fi
    
    echo "NFS server setup complete"
}

# Function to setup NFS client
setup_nfs_client() {
    local nfs_server=$1
    
    echo "Setting up NFS client for server: $nfs_server"
    
    # Install dependencies
    install_common_deps
    
    # Install rclone
    curl https://rclone.org/install.sh | sudo bash
    
    # Download and install scratch-dna
    wget https://github.com/FredHutch/sc-benchmark/raw/master/bin/scratch-dna-go_linux-amd64 -O scratch-dna
    chmod +x scratch-dna
    sudo mv scratch-dna /usr/local/bin/
    
    # Create mount point
    sudo mkdir -p /mnt/test
    
    # Mount NFS share
    sudo mount -t nfs4  -o noac,sync,actimeo=0 ${nfs_server}:/opt /mnt/test
    
    # Run benchmark
    echo "Running scratch-dna benchmark, writing 1MB files ..."
    scratch-dna -v -p 1 4096 1048576 1 /mnt/test/temp/
}

# Main script
if [ $# -eq 0 ]; then
    setup_nfs_server
else
    setup_nfs_client "$1"
fi

