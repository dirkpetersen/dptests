#! /bin/bash

CEPH_RELEASE=19.2.2  # replace this with the active release

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
#/usr/local/bin/cephadm bootstrap --allow-fqdn-hostname --mon-ip ${INTERNAL_IP} --ssh-user rocky
/usr/local/bin/cephadm bootstrap --mon-ip ${INTERNAL_IP} --ssh-user rocky
