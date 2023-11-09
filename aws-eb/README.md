# AWS-EB (Easybuild in AWS)

This tool builds HPC software using the Easybuild (EB) Framework. It attempts to compile the latest version of all packages and then tar and upload these builds to AWS S3. The EB root (EASYBUILD_PREFIX) is currently set to /opt/eb.

## Try it

Make sure you have your AWS creds in ~/.aws/credentials and your AWS region is set. First run the config and then launch --list to see what cpu types are there. 

```
./aws-eb.py config
./aws-eb.py launch --list
./aws-eb.py launch --cpu-type epyc-gen-4 --bio-only
```

## Build Workflow 

This is happening behind the scenes: 

1. Launch AWS instance 
1. Attach 1TB EBS volumne to /opt/eb
1. Upload and launch bootstrap.sh script and other configs
1. Install basic software (Python3, Apptainer, etc) 
1. Launch aws-eb script with `--build` option and same parameters as launched first
1. Download all tarred binaries from S3 and unpack them and download sources and modules
1. Loop through all *.eb files (except `__archive__`) and for each eb file:
   1. check for allowed toolchains and see if --bio-only
   1. install all osdependencies via dnf or apt
   1. download software from shared S3 bucket 
   1. set all files under sources/generic to executable 
   1. untar .eb.tar.gz files 
   1. build the latest version of each software using `eb --robot`
   1. tar new software to .eb.tar.gz files 
   1. upload software to shared S3 bucket 
1. run `aws-eb.py ssh --terminate <instance>` after the all easyconfigs are processed.

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

The most cost efficient instance type is not clear yet. I started with c7a.xlarge with 4 vcpus and 8GB RAM. It may not make sense to use a larger instance type as there are long periods of time where only a single vcpu is running. At the tail end it installs R packages for hours which is limited to a single vcpu. Perhaps run just more instances  

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

## allowed toolchains 

These are currently the only allowed toolchains 

`allowed_toolchains = ['system', 'GCC', 'GCCcore', 'foss', 'fosscuda']`

## CLI

```
usage: aws-eb  [-h] [--debug] [--profile AWSPROFILE] [--version]
               {config,cnf,launch,lau,download,dld,ssh,scp} ...

A (mostly) for building stuff on AWS after finding folders in the file system that are worth
archiving.

positional arguments:
  {config,cnf,launch,lau,download,dld,ssh,scp}
                        sub-command help
    config (cnf)        Bootstrap the configuration, install dependencies and setup your
                        environment. You will need to answer a few questions about your cloud and
                        hpc setup.
    launch (lau)        Launch EC2 instance, build new Easybuild packages and upload them to S3
    download (dld)      Download built eb packages to /opt/eb
    ssh (scp)           Login to an AWS EC2 instance 

optional arguments:
  -h, --help            show this help message and exit
  --debug, -d           verbose output for all commands
  --profile AWSPROFILE, -p AWSPROFILE
                        which AWS profile in ~/.aws/ should be used. default="aws"
  --version, -v         print AWS-EB and Python version info
```


```
./aws-eb.py config --help
Error: EasyBuild not found. Please install it first.
usage: aws-eb config [-h] [--monitor <email@address.org>]

optional arguments:
  -h, --help            show this help message and exit
  --monitor <email@address.org>, -m <email@address.org>
                        setup aws-eb as a monitoring cronjob on an ec2 instance and notify an email address
```

GPU is not yet implemented 

```
./aws-eb.py launch --help
usage: aws-eb launch [-h] [--instance-type INSTANCETYPE] [--gpu-type GPUTYPE] [--cpu-type CPUTYPE]
                     [--list] [--vcpus VCPUS] [--mem MEM] [--monitor] [--build] [--bio-only]

optional arguments:
  -h, --help            show this help message and exit
  --instance-type INSTANCETYPE, -i INSTANCETYPE
                        The EC2 instance type is auto-selected, but you can pick any other type here
  --gpu-type GPUTYPE, -g GPUTYPE
                        run --list to see available GPU types
  --cpu-type CPUTYPE, -c CPUTYPE
                        run --list to see available CPU types
  --list, -l            List CPU and GPU types
  --vcpus VCPUS, -v VCPUS
                        Number of cores to be allocated for the machine. (default=4)
  --mem MEM, -m MEM     GB Memory allocated to instance  (default=8)
  --monitor, -n         Monitor EC2 server for cost and idle time.
  --build, -b           Build the Easybuild packages on current system.
  --bio-only, -o        Build only life sciences (bio) packages and dependencies.
```

just run `./aws-eb.py ssh` to connect to the EC2 instance you created last

```
./aws-eb.py ssh --help
usage: aws-eb ssh [-h] [--list] [--terminate <hostname>] [sshargs ...]

positional arguments:
  sshargs               multiple arguments to ssh/scp such as hostname or user@hostname oder folder

optional arguments:
  -h, --help            show this help message and exit
  --list, -l            List running AWS-EB EC2 instances
  --terminate <hostname>, -t <hostname>
                        Terminate EC2 instance with this public IP Address.
```


Standalone Download is not yet implemented 

```
./aws-eb.py download --help
usage: aws-eb download [-h] [--target <target_folder>]

optional arguments:
  -h, --help            show this help message and exit
  --target <target_folder>, -t <target_folder>
                        Download to other folder than default
```

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
