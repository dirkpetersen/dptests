# AWS-EB (Easybuild in AWS)

This tool builds HPC Software using the Easybuild (EB) Framework. It attempts to build the latest version of all packages and then tar and upload these builds to AWS S3. The EB root (EASYBUILD_PREFIX) is currently set to /opt/eb.

## Build Workflow 

1. Launch AWS instance 
1. Upload and launch bootstrap.sh script 
1. Install basic software (Python3, Apptainer, etc) 
1. Launch aws-eb script with --build option and same parameters as launched first
1. Download all tarred binaries from S3 and unpack them and download sources and modules
1. Loop through all *.eb files (except __archive__) and build the latest version of each software
1. Tar up results after each build and upload to S3 
1. Downloads happen before AND after each eb --robot execution. This allows multiple build machines to be running in parallel building for the same architecture while they are synced using S3

## S3 Folder Layout 

Tarred binaries and modules are copied to a platform specific folder (e.g amzn-2023_epyc-gen-4/software) and sources are copies to a shared folder that all platforms use.

![image](https://github.com/dirkpetersen/dptests/assets/1427719/91b2542e-81c4-4859-ad9e-d6ed7b7ac1c3)


## Instance Mapping 

Each CPU family or GPU type is mapped to all AWS instance families that have this CPU family or GPU installed.
This will allow to pick any spot instance that has a certain compatible hardware configuration. 

```
        self.cpu_types = {
            "graviton-2": ['m6g','c6g', 'c6gn', 't4g' ,'g5g'],
            "graviton-3": ['m7g', 'c7g', 'c7gn'],
            "epyc-gen-1": ['t3a'],
            "epyc-gen-2": ['c5a', 'm5a', 'r5a', 'g4ad', 'p4', 'inf2', 'g5'],
            "epyc-gen-3": ['c6a', 'm6a', 'r6a', 'p5'],
            "epyc-gen-4": ['c7a', 'm7a', 'r7a'],
            "xeon-gen-1": ['c4', 'm4', 't2', 'r4', 'p3' ,'p2', 'f1', 'g3'],
            "xeon-gen-2": ['c5', 'm5', 'c5n', 'm5n', 'm5zn', 'r5', 't3', 't3n', 'dl1', 'inf1', 'g4dn', 'vt1'],
            "xeon-gen-3": ['c6i', 'm6i', 'm6in', 'c6in', 'r6i', 'r6id', 'r6idn', 'r6in', 'trn1'],
            "xeon-gen-4": ['c7i', 'm7i', 'm7i-flex'],
            "core-i7-mac": ['mac1']
        }
        self.gpu_types = {
            "h100": 'p5',
            "a100": 'p4',
            "v100": 'p3',  
            "k80": 'p2',
            "gaudi": 'dl1',
            "trainium": 'trn1',
            "inferentia2": 'inf2',
            "inferentia1": 'inf1',
            "t4g": 'g5g',
            "a10g": 'g5',
            "t4": 'g4dn',
            "v520": 'g4ad',
            "m60": 'g3',
            "fpga": 'f1',
            "u30": 'vt1'            
        }
```

## build-machine 

The most cost efficient instance type is not clear yet. I started with c7a.xlarge with 4 cpus and 8GB RAM. It is not clear if it would make sense to use a larger instance type As there are long periods of time where only a single vcpu is running. At the tail end it installs R packages for hours which is limited to a single CPU 

![image](https://github.com/dirkpetersen/dptests/assets/1427719/973f973d-14f0-41af-8f7e-838d59ad58f4)



It has this envionment. 

```
(base) ec2-user@froster:/opt/eb$ cat ~/.easybuildrc
#! /bin/bash

source /usr/local/lmod/lmod/init/bash
export MODULEPATH=/opt/eb/modules/all
#
export EASYBUILD_JOB_CORES=4
export EASYBUILD_CUDA_COMPUTE_CAPABILITIES=7.5,8.0,8.6,9.0
# export EASYBUILD_BUILDPATH=/dev/shm/$USER
export EASYBUILD_PREFIX=/opt/eb
export EASYBUILD_JOB_OUTPUT_DIR=$EASYBUILD_PREFIX/batch-output
export EASYBUILD_JOB_BACKEND=Slurm
export EASYBUILD_PARALLEL=16
# export EASYBUILD_GITHUB_USER=${WHOAMI}
export  EASYBUILD_UPDATE_MODULES_TOOL_CACHE=True
export  EASYBUILD_ROBOT_PATHS=/home/ec2-user/.local/easybuild/easyconfigs:/opt/eb/fh/fh_easyconfigs
```


## eb toolchains 

allowed_toolchains = ['GCC', 'GCCcore', 'SYSTEM', 'foss', 'fosscuda']

/home/ec2-user/miniconda3/easybuild/easyconfigs/x/xESMF/xESMF-0.3.0-intel-2020b.eb:toolchain = {'name': 'intel', 'version': '2020b'}
/home/ec2-user/miniconda3/easybuild/easyconfigs/x/xarray/xarray-0.16.2-fosscuda-2020b.eb:toolchain = {'name': 'fosscuda', 'version': '2020b'}
/home/ec2-user/miniconda3/easybuild/easyconfigs/x/xarray/xarray-2022.6.0-foss-2022a.eb:toolchain = {'name': 'foss', 'version': '2022a'}
/home/ec2-user/miniconda3/easybuild/easyconfigs/x/xarray/xarray-2023.4.2-gfbf-2022b.eb:toolchain = {'name': 'gfbf', 'version': '2022b'}
/home/ec2-user/miniconda3/easybuild/easyconfigs/x/xbitmaps/xbitmaps-1.1.1.eb:toolchain = SYSTEM
/home/ec2-user/miniconda3/easybuild/easyconfigs/x/xclip/xclip-0.13-GCCcore-11.3.0.eb:toolchain = {'name': 'GCCcore', 'version': '11.3.0'}
/home/ec2-user/miniconda3/easybuild/easyconfigs/x/xmlf90/xmlf90-1.5.4-GCC-10.2.0.eb:toolchain = {'name': 'GCC', 'version': '10.2.0'}
/home/ec2-user/miniconda3/easybuild/easyconfigs/x/xmlf90/xmlf90-1.5.4-iccifort-2019.5.281.eb:toolchain = {'name': 'iccifort', 'version': '2019.5.281'}



## supported OS 

```
(base) ec2-user@froster:~$ cat /etc/os-release
NAME="Amazon Linux"
VERSION="2023"
ID="amzn"
ID_LIKE="fedora"
VERSION_ID="2023"
PLATFORM_ID="platform:al2023"
PRETTY_NAME="Amazon Linux 2023"
ANSI_COLOR="0;33"
CPE_NAME="cpe:2.3:o:amazon:amazon_linux:2023"
HOME_URL="https://aws.amazon.com/linux/"
BUG_REPORT_URL="https://github.com/amazonlinux/amazon-linux-2023"
SUPPORT_END="2028-03-01"
```

```
petersen@rhino03:~$ cat /etc/os-release
NAME="Ubuntu"
VERSION="18.04.6 LTS (Bionic Beaver)"
ID=ubuntu
ID_LIKE=debian
PRETTY_NAME="Ubuntu 18.04.6 LTS"
VERSION_ID="18.04"
HOME_URL="https://www.ubuntu.com/"
SUPPORT_URL="https://help.ubuntu.com/"
BUG_REPORT_URL="https://bugs.launchpad.net/ubuntu/"
PRIVACY_POLICY_URL="https://www.ubuntu.com/legal/terms-and-policies/privacy-policy"
VERSION_CODENAME=bionic
UBUNTU_CODENAME=bionic
```

```
[peterdir@exanode-08-1 ~]$ cat /etc/os-release
NAME="CentOS Linux"
VERSION="7 (Core)"
ID="centos"
ID_LIKE="rhel fedora"
VERSION_ID="7"
PRETTY_NAME="CentOS Linux 7 (Core)"
ANSI_COLOR="0;31"
CPE_NAME="cpe:/o:centos:centos:7"
HOME_URL="https://www.centos.org/"
BUG_REPORT_URL="https://bugs.centos.org/"
CENTOS_MANTISBT_PROJECT="CentOS-7"
CENTOS_MANTISBT_PROJECT_VERSION="7"
REDHAT_SUPPORT_PRODUCT="centos"
REDHAT_SUPPORT_PRODUCT_VERSION="7"
```

```
dp@grammy:~/gh/dptests$ cat /etc/os-release
PRETTY_NAME="Debian GNU/Linux 11 (bullseye)"
NAME="Debian GNU/Linux"
VERSION_ID="11"
VERSION="11 (bullseye)"
VERSION_CODENAME=bullseye
ID=debian
HOME_URL="https://www.debian.org/"
SUPPORT_URL="https://www.debian.org/support"
BUG_REPORT_URL="https://bugs.debian.org/"
```
