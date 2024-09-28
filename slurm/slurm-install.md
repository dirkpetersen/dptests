# How to install Slurm on a single machine 

This will install the Slurm Workloadmanager on a Rocky9 (RHEL9) Machine as localhost :

```
dnf install epel-release
dnf config-manager --set-enabled crb
dnf update
dnf groupinstall "Development Tools"
dnf install libjwt-devel hwloc-devel bpftool dbus-devel hdf5-devel man2html curl-devel json-c-devel json-devel libibmad-devel libibumad-devel libvpl-devel freeipmi-devel lua-devel munge-devel mariadb-devel numactl-devel pam-devel pmix-devel readline-devel http-parser-devel gtk2-devel http-parser-devel libyaml-devel perl-ExtUtils-MakeMaker

wget https://download.schedmd.com/slurm/slurm-24......tar.bz2
rpmbuild -ta slurm-24.....bz2 --with slurmrestd

rpm -i /root/rpmbuild/RPMS/x86_64/slurm-24.....rpm
rpm -i /root/rpmbuild/RPMS/x86_64/slurm-slurmctld-24.....rpm
rpm -i /root/rpmbuild/RPMS/x86_64/slurm-slurmd-24.....rpm
rpm -i /root/rpmbuild/RPMS/x86_64/slurm-devel-24.....rpm
rpm -i /root/rpmbuild/RPMS/x86_64/slurm-perlapi-24.....rpm
rpm -i /root/rpmbuild/RPMS/x86_64/slurm-example-configs-24.....rpm

sed -e 's/SlurmctldHost=linux0/SlurmctldHost=localhost/' -e 's/NodeName=linux\[1-32\]/NodeName=localhost/' /etc/slurm/slurm.conf.example > /etc/slurm/slurm.conf
cp /etc/slurm/cgroup.conf.example /etc/slurm/cgroup.conf

useradd slurm
mkdir /var/spool/slurmd
mkdir /var/spool/slurmctld
chown slurm /var/spool/slurmd
chown slurm /var/spool/slurmctld

/usr/sbin/create-munge-key

systemctl enable munge slurmctld slurmd
systemctl restart munge slurmctld slurmd

```

now you can test something like

```
srun hostname
```




