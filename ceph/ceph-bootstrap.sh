CEPH_RELEASE=19.2.1 # replace this with the active release
curl -o /usr/local/bin/cephadm https://download.ceph.com/rpm-${CEPH_RELEASE}/el9/noarch/cephadm
chmod +x /usr/local/bin/cephadm
cephadm bootstrap --allow-fqdn-hostname --mon-ip ${MON_IP}
