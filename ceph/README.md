# Ceph Cluster Deployment Toolkit

This repository contains a comprehensive set of shell scripts for automating the deployment and management of Ceph storage clusters on AWS EC2 infrastructure.

## Overview

The toolkit provides end-to-end automation for:
- AWS EC2 instance provisioning with proper security groups
- EBS volume management and attachment
- Ceph cluster bootstrapping and node orchestration
- Cross-host SSH key distribution and file management
- Lambda-based monitoring and cleanup functions

## Core Scripts

### 🚀 Main Deployment Scripts

#### `ec2-create-instances.sh`
**Primary cluster deployment orchestrator**

Creates and configures AWS EC2 instances for Ceph cluster deployment.

**Key Features:**
- **Instance Management**: Discovers existing instances, launches missing ones
- **Storage Configuration**: Attaches configurable EBS volumes (currently 7x 125GB ST1)
- **Network Setup**: Configures security groups with automatic SSH self-access rules
- **DNS Integration**: Registers Route 53 A records for cluster nodes
- **Ceph Orchestration**: Bootstraps first node, then adds additional nodes to cluster

**Usage:**
```bash
# Launch 3-node cluster
./ec2-create-instances.sh 3

# Launch single node  
./ec2-create-instances.sh 1
```

**Configuration (Environment Variables):**
- `AWS_REGION="us-west-2"` - AWS region
- `EC2_TYPE="c8gd.large"` - Instance type (c8gd.large/c5ad.large/c7gd.medium)
- `AMI_IMAGE="ami-03be04a3da3a40226"` - Rocky Linux 9 ARM64 AMI
- `EC2_SECURITY_GROUPS="SSH-HTTP-ICMP ceph-cluster-sg"` - Security groups
- `EBS_QTY="7"` - Number of EBS volumes per instance
- `EBS_SIZE="125"` - Size of each EBS volume (GB)
- `EBS_TYPE="st1"` - EBS volume type (st1/gp3/io2)

#### `bootstrap-node.sh`
**Ceph cluster bootstrap and node management**

Handles Ceph cluster initialization and node integration with intelligent device detection.

**Operational Modes:**
- **`first`** - Bootstrap new Ceph cluster (MON node)
- **`others`** - Prepare additional nodes for cluster joining
- **`join_cluster`** - Add nodes to existing cluster via orchestrator
- **`create_osds`** - Create Object Storage Daemons on target hosts

**Key Features:**
- **Smart Device Detection**: Retry mechanism with intelligent analysis
  - Distinguishes between available, in-use, and rejected devices
  - Stops retrying permanently rejected devices (partition tables, size issues)
  - Continues retrying for legitimate timing issues on new nodes
- **Comprehensive Status Reporting**: Detailed device analysis and error messages
- **HDD/SSD Optimization**: Configurable ratios for shared DB/WAL on SSDs
- **Cluster Integration**: Waits for proper orchestrator connectivity

**Usage:**
```bash
# Bootstrap first node
sudo bash bootstrap-node.sh first

# Prepare additional node
sudo bash bootstrap-node.sh others

# Add node to cluster (run from first node)
sudo TARGET_HOSTNAME=ceph-test-2 TARGET_INTERNAL_IP=10.0.1.100 bash bootstrap-node.sh join_cluster

# Create OSDs (run from first node)  
sudo TARGET_HOSTNAME=ceph-test-2 bash bootstrap-node.sh create_osds
```

### 🔧 Infrastructure Support Scripts

#### `ec2-create-ceph-security-group.sh`
**AWS security group creation for Ceph clusters**

Creates and configures security groups with all necessary Ceph ports and protocols.

**Features:**
- Creates `ceph-cluster-sg` security group
- Configures ingress rules for Ceph services:
  - MON: 3300, 6789
  - MGR: 9283 (dashboard), 8765, 8443
  - OSD: 6800-7300
  - MDS: 6800
  - RGW: 7480, 8080
  - SSH: 22 (from security group itself for orchestration)

#### `ebs-create-attach.sh`
**EBS volume provisioning and attachment**

Creates and attaches EBS volumes to EC2 instances for Ceph storage.

**Usage:**
```bash
# Attach 6x 125GB ST1 volumes to instance
./ebs-create-attach.sh i-1234567890abcdef0 st1 125 6
```

**Parameters:**
1. Instance ID
2. Volume type (st1, gp3, io2)
3. Volume size (GB)
4. Number of volumes

#### `ebs-delete-unused.sh`
**EBS volume cleanup utility**

Identifies and deletes unattached EBS volumes to prevent cost accumulation.

**Features:**
- Lists all unattached volumes
- Provides deletion commands
- Safety checks to prevent accidental deletion of attached volumes

### 🔄 Utility Scripts

#### `remote-file-copy.sh`
**Cross-host file transfer and SSH key distribution**

Provides functions for secure file transfers between remote hosts using the local machine as a relay.

**Key Functions:**
- **`remote_file_copy()`** - Copy files between hosts preserving permissions
- **`remote_ssh_key_copy()`** - Distribute SSH keys to authorized_keys files

**Features:**
- Works through local machine (no direct host-to-host connections needed)
- Preserves file permissions and ownership
- Handles sudo access and protected directories
- Creates target directories automatically

**Usage:**
```bash
# Copy configuration file
remote_file_copy 'user@host1:/etc/ceph/ceph.conf' 'user@host2:/etc/ceph/ceph.conf'

# Copy SSH key
remote_ssh_key_copy 'user@host1:/etc/ceph/ceph.pub' 'user@host2' 'root'
```

#### `deploy-lambda.sh`
**AWS Lambda deployment automation**

Packages and deploys Lambda functions with proper IAM roles and permissions.

**Features:**
- Creates deployment packages
- Configures IAM roles and policies
- Handles Lambda function updates
- Sets up CloudWatch event triggers

## Configuration Files

### `ec2-cloud-init.txt`
Cloud-init configuration for EC2 instances. Handles:
- Package installation (podman, lvm2, etc.)
- System configuration
- Initial setup tasks

## Deployment Workflow

1. **Security Setup**: Run `ec2-create-ceph-security-group.sh`
2. **Cluster Launch**: Run `ec2-create-instances.sh [num_instances]`
3. **Automatic Process**:
   - Launches EC2 instances
   - Attaches EBS volumes
   - Bootstraps first node with Ceph
   - Adds additional nodes to cluster
   - Creates OSDs with optimal HDD/SSD ratios

## Architecture

- **First Node**: Acts as cluster orchestrator (runs MON, MGR services)
- **Additional Nodes**: Join cluster and provide OSD services
- **Storage**: ST1 volumes for cost-effective bulk storage
- **Network**: Private cluster communication + public management access

## Current Configuration

- **Instance Type**: c8gd.large (ARM64, NVMe SSD + network-optimized)
- **AMI**: Rocky Linux 9 ARM64 (`ami-03be04a3da3a40226`)
- **Storage**: 7x 125GB ST1 volumes per node
- **Ceph Version**: 19.2.2 (latest stable)

## Monitoring & Maintenance

- **Cluster Status**: `sudo /usr/local/bin/cephadm shell -- ceph -s`
- **Device Status**: `sudo /usr/local/bin/cephadm shell -- ceph orch device ls`
- **OSD Status**: `sudo /usr/local/bin/cephadm shell -- ceph osd tree`
- **Clean Unused Volumes**: Run `ebs-delete-unused.sh` periodically

## Troubleshooting

### Mount Ceph Container
```bash
sudo /usr/local/bin/cephadm shell
```

### Manual OSD Creation
If the Web UI has issues with AWS device detection:

```bash
# Single device
ceph orch daemon add osd ceph-test-1:/dev/nvme2n1

# All available devices
ceph orch daemon add osd <hostname> --all-available-devices
```

### Device Detection Issues
The toolkit includes intelligent device detection that handles:
- **Timing Issues**: Newly added nodes need time for orchestrator inventory
- **Rejected Devices**: Devices with partitions, wrong size, or hardware issues
- **In-Use Devices**: Already configured devices (normal condition)

Check device status with detailed analysis:
```bash
ceph orch device ls --format json
```

## Security Features

- **Network Isolation**: Security groups restrict access to necessary ports only
- **SSH Key Management**: Automated distribution of Ceph orchestrator keys
- **Least Privilege**: IAM roles with minimal required permissions
- **Encrypted Storage**: EBS volumes support encryption at rest

## Cost Optimization

- **ST1 Volumes**: Throughput-optimized for bulk storage workloads
- **Instance Right-sizing**: c8gd.large provides optimal CPU/memory/network balance
- **Volume Cleanup**: Automated detection and cleanup of unused EBS volumes
- **Spot Instance Support**: Can be configured for non-production workloads