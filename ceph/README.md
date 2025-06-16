
# Troubleshooting 

## mount ceph container 

- sudo /usr/local/bin/cephadm shell

## Create OSD, no devices were found 

Web UI issue with AWS, use instead 

```
$ ceph orch daemon add osd ceph-test-1:/dev/nvme2n1
Created osd(s) 0 on host 'ceph-test-1'
$ ceph orch daemon add osd ceph-test-1:/dev/nvme3n1
Created osd(s) 1 on host 'ceph-test-1'
```

or 

```
ceph orch daemon add osd <hostname> --all-available-devices
```