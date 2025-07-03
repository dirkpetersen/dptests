#!/bin/bash

# Variables - Modify these as needed
VPC_ID="vpc-xxxxxxxxx"  # Replace with your VPC ID
SG_NAME="ceph-cluster-sg"
SG_DESCRIPTION="Security Group for Ceph Cluster Inter-node Communication"
REGION="us-west-2"  # Replace with your region

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${YELLOW}Creating Ceph Cluster Security Group...${NC}"

# Create the security group
SG_ID=$(aws ec2 create-security-group \
    --group-name "$SG_NAME" \
    --description "$SG_DESCRIPTION" \
    --vpc-id "$VPC_ID" \
    --region "$REGION" \
    --query 'GroupId' \
    --output text)

if [ $? -eq 0 ]; then
    echo -e "${GREEN}Security Group created successfully: $SG_ID${NC}"
else
    echo -e "${RED}Failed to create security group${NC}"
    exit 1
fi

# Add tags to the security group
aws ec2 create-tags \
    --resources "$SG_ID" \
    --tags Key=Name,Value="$SG_NAME" Key=Purpose,Value="Ceph Cluster" \
    --region "$REGION"

echo -e "${YELLOW}Adding ingress rules...${NC}"

# Function to add ingress rule
add_ingress_rule() {
    local port=$1
    local protocol=$2
    local description=$3
    
    aws ec2 authorize-security-group-ingress \
        --group-id "$SG_ID" \
        --protocol "$protocol" \
        --port "$port" \
        --source-group "$SG_ID" \
        --region "$REGION" \
        --rule-description "$description" >/dev/null 2>&1
    
    if [ $? -eq 0 ]; then
        echo -e "${GREEN}✓ Added rule: $description (Port: $port)${NC}"
    else
        echo -e "${RED}✗ Failed to add rule: $description (Port: $port)${NC}"
    fi
}

# Function to add port range rule
add_port_range_rule() {
    local from_port=$1
    local to_port=$2
    local protocol=$3
    local description=$4
    
    aws ec2 authorize-security-group-ingress \
        --group-id "$SG_ID" \
        --protocol "$protocol" \
        --port "$from_port-$to_port" \
        --source-group "$SG_ID" \
        --region "$REGION" \
        --rule-description "$description" >/dev/null 2>&1
    
    if [ $? -eq 0 ]; then
        echo -e "${GREEN}✓ Added rule: $description (Ports: $from_port-$to_port)${NC}"
    else
        echo -e "${RED}✗ Failed to add rule: $description (Ports: $from_port-$to_port)${NC}"
    fi
}

# Ceph Monitor ports
add_ingress_rule "3300" "tcp" "Ceph Monitor (new default)"
add_ingress_rule "6789" "tcp" "Ceph Monitor (legacy)"

# OSD/MGR/MDS port range
add_port_range_rule "6800" "7300" "tcp" "Ceph OSD/MGR/MDS services"

# Ceph Dashboard ports
add_ingress_rule "8080" "tcp" "Ceph Dashboard HTTP"
add_ingress_rule "8443" "tcp" "Ceph Dashboard HTTPS"

# RGW ports
add_ingress_rule "7480" "tcp" "Ceph RGW HTTP"
add_ingress_rule "7481" "tcp" "Ceph RGW HTTPS"
add_ingress_rule "80" "tcp" "HTTP (RGW alternative)"
add_ingress_rule "443" "tcp" "HTTPS (RGW alternative)"

# iSCSI Gateway (if used)
add_ingress_rule "3260" "tcp" "iSCSI target"

# SSH for management (optional but recommended)
add_ingress_rule "22" "tcp" "SSH management"

# NTP for time synchronization (important for Ceph)
add_ingress_rule "123" "udp" "NTP time sync"

# Add ICMP for ping/troubleshooting
aws ec2 authorize-security-group-ingress \
    --group-id "$SG_ID" \
    --protocol "icmp" \
    --port "-1" \
    --source-group "$SG_ID" \
    --region "$REGION" \
    --rule-description "ICMP for ping/troubleshooting" >/dev/null 2>&1

if [ $? -eq 0 ]; then
    echo -e "${GREEN}✓ Added rule: ICMP for ping/troubleshooting${NC}"
else
    echo -e "${RED}✗ Failed to add rule: ICMP${NC}"
fi

echo -e "${YELLOW}Summary:${NC}"
echo -e "Security Group ID: ${GREEN}$SG_ID${NC}"
echo -e "Security Group Name: ${GREEN}$SG_NAME${NC}"
echo -e "VPC ID: ${GREEN}$VPC_ID${NC}"
echo -e "Region: ${GREEN}$REGION${NC}"

echo -e "\n${YELLOW}To attach this security group to your Ceph nodes, use:${NC}"
echo -e "${GREEN}aws ec2 modify-instance-attribute --instance-id <instance-id> --groups $SG_ID --region $REGION${NC}"

echo -e "\n${YELLOW}To view the security group rules:${NC}"
echo -e "${GREEN}aws ec2 describe-security-groups --group-ids $SG_ID --region $REGION${NC}"

echo -e "\n${GREEN}Ceph cluster security group setup completed!${NC}"