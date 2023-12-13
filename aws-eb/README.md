# AWS-EB (Easybuild in AWS)

This tool builds HPC software using the Easybuild (EB) Framework. It attempts to compile the latest version of all packages and then tar and upload these builds to AWS S3. The EB root (EASYBUILD_PREFIX) is currently set to /opt/eb.

## Try it

Make sure you have your AWS creds in ~/.aws/credentials and your AWS region is set. First run the config and then launch --list to see what cpu types are there. 

```
./aws-eb.py config
./aws-eb.py launch --list
./aws-eb.py launch --cpu-type epyc-gen-4 --include bio,lib
```

## Build Workflow 

This is happening behind the scenes: 

1. Launch AWS instance 
1. Attach 750 GB EBS volumne to /opt/eb
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


## after 1 day 

after 24h of building I see this: 

```
ec2-user@aws-eb:~$ ml avail

--------------------------------------- /opt/eb/modules/all ----------------------------------------
   ADMIXTURE/1.3.0                                     XML-LibXML/2.0207-GCCcore-11.3.0
   ATK/2.38.0-GCCcore-11.3.0                           XPLOR-NIH/3.4-Linux_x86_64
   Archive-Zip/1.68-GCCcore-11.3.0                     XZ/5.2.5-GCCcore-11.2.0
   Arrow/8.0.0-foss-2022a                              XZ/5.2.5-GCCcore-11.3.0             (D)
   Autoconf/2.71-GCCcore-11.2.0                        Xvfb/21.1.3-GCCcore-11.3.0
   Autoconf/2.71-GCCcore-11.3.0                        ZeroMQ/4.3.4-GCCcore-11.3.0
   Autoconf/2.71                               (D)     Zip/3.0-GCCcore-11.3.0
   Automake/1.16.4-GCCcore-11.2.0                      arrow-R/8.0.0-foss-2022a-R-4.2.1
   Automake/1.16.5-GCCcore-11.3.0                      binutils/2.32
   Automake/1.16.5                             (D)     binutils/2.34
   Autotools/20210726-GCCcore-11.2.0                   binutils/2.35
   Autotools/20220317-GCCcore-11.3.0                   binutils/2.36.1
   Autotools/20220317                          (D)     binutils/2.37-GCCcore-11.2.0
   BLAST/2.2.26-Linux_x86_64                           binutils/2.37
   BLIS/0.9.0-GCC-11.3.0                               binutils/2.38-GCCcore-11.3.0
   BeautifulSoup/4.10.0-GCCcore-11.3.0                 binutils/2.38
   Bio-DB-HTS/3.01-GCC-11.3.0                          binutils/2.39-GCCcore-12.2.0
   BioPerl/1.7.8-GCCcore-11.3.0                        binutils/2.39                       (D)
   Biopython/1.79-foss-2022a                           bzip2/1.0.8-GCCcore-11.2.0
   Bison/3.3.2                                         bzip2/1.0.8-GCCcore-11.3.0          (D)
   Bison/3.7.6-GCCcore-11.2.0                          cURL/7.78.0-GCCcore-11.2.0
   Bison/3.8.2-GCCcore-11.3.0                          cURL/7.83.0-GCCcore-11.3.0          (D)
   Bison/3.8.2-GCCcore-12.2.0                          cairo/1.17.4-GCCcore-11.3.0
   Bison/3.8.2                                 (D)     expat/2.4.1-GCCcore-11.2.0
   Blosc/1.21.3-GCCcore-11.3.0                         expat/2.4.8-GCCcore-11.3.0
   Blosc2/2.6.1-GCCcore-11.3.0                         expat/2.4.9-GCCcore-12.2.0          (D)
   Boost/1.79.0-GCC-11.3.0                             flex/2.6.4-GCCcore-11.2.0
   Brotli/1.0.9-GCCcore-11.2.0                         flex/2.6.4-GCCcore-11.3.0
   Brotli/1.0.9-GCCcore-11.3.0                 (D)     flex/2.6.4-GCCcore-12.2.0
   CMake/3.21.1-GCCcore-11.2.0                         flex/2.6.4                          (D)
   CMake/3.23.1-GCCcore-11.3.0                         fontconfig/2.14.0-GCCcore-11.3.0
   CMake/3.24.3-GCCcore-11.3.0                 (D)     foss/2022a
   Compress-Raw-Zlib/2.202-GCCcore-11.3.0              freetype/2.11.0-GCCcore-11.2.0
   DB/18.1.40-GCCcore-11.2.0                           freetype/2.12.1-GCCcore-11.3.0      (D)
   DB/18.1.40-GCCcore-11.3.0                           gettext/0.21-GCCcore-11.2.0
   DB/18.1.40-GCCcore-12.2.0                   (D)     gettext/0.21-GCCcore-11.3.0
   DB_File/1.858-GCCcore-11.3.0                        gettext/0.21                        (D)
   Doxygen/1.9.4-GCCcore-11.3.0                        git/2.33.1-GCCcore-11.2.0-nodocs
   Eigen/3.4.0-GCCcore-11.3.0                          git/2.36.0-GCCcore-11.3.0-nodocs    (D)
   FFTW.MPI/3.3.10-gompi-2022a                         gompi/2022a
   FFTW/3.3.10-GCC-11.3.0                              googletest/1.11.0-GCCcore-11.3.0
   FLAC/1.3.4-GCCcore-11.3.0                           gperf/3.1-GCCcore-11.2.0
   Flask/2.2.2-GCCcore-11.3.0                          gperf/3.1-GCCcore-11.3.0            (D)
   FlexiBLAS/3.2.0-GCC-11.3.0                          groff/1.22.4-GCCcore-11.2.0
   FriBidi/1.0.12-GCCcore-11.3.0                       groff/1.22.4-GCCcore-11.3.0
   GCC/11.2.0                                          groff/1.22.4-GCCcore-12.2.0         (D)
   GCC/11.3.0                                          gzip/1.12-GCCcore-11.3.0
   GCC/12.2.0                                  (D)     help2man/1.48.3-GCCcore-11.2.0
   GCCcore/11.2.0                                      help2man/1.49.2-GCCcore-11.3.0
   GCCcore/11.3.0                                      help2man/1.49.2-GCCcore-12.2.0      (D)
   GCCcore/12.2.0                              (D)     hwloc/2.7.1-GCCcore-11.3.0
   GDAL/3.5.0-foss-2022a                               hypothesis/6.46.7-GCCcore-11.3.0
   GEOS/3.10.3-GCC-11.3.0                              intltool/0.51.0-GCCcore-11.2.0
   GLPK/5.0-GCCcore-11.3.0                             intltool/0.51.0-GCCcore-11.3.0      (D)
   GLib/2.72.1-GCCcore-11.3.0                          jbigkit/2.1-GCCcore-11.3.0
   GMP/6.2.1-GCCcore-11.2.0                            jemalloc/5.3.0-GCCcore-11.3.0
   GMP/6.2.1-GCCcore-11.3.0                    (D)     libGLU/9.0.2-GCCcore-11.3.0
   GObject-Introspection/1.72.0-GCCcore-11.3.0         libaec/1.0.6-GCCcore-11.3.0
   GSL/2.7-GCC-11.2.0                                  libaio/0.3.112-GCCcore-11.3.0
   GSL/2.7-GCC-11.3.0                          (D)     libarchive/3.5.1-GCCcore-11.2.0
   GTK2/2.24.33-GCCcore-11.3.0                         libarchive/3.6.1-GCCcore-11.3.0     (D)
   Gdk-Pixbuf/2.42.8-GCCcore-11.3.0                    libdeflate/1.10-GCCcore-11.3.0
   Ghostscript/9.56.1-GCCcore-11.3.0                   libdrm/2.4.110-GCCcore-11.3.0
   GitPython/3.1.27-GCCcore-11.3.0                     libevent/2.1.12-GCCcore-11.3.0
   HDF/4.2.15-GCCcore-11.3.0                           libfabric/1.15.1-GCCcore-11.3.0
   HDF5/1.12.2-gompi-2022a                             libffi/3.4.2-GCCcore-11.2.0
   HTSlib/1.15.1-GCC-11.3.0                            libffi/3.4.2-GCCcore-11.3.0         (D)
   HarfBuzz/4.2.1-GCCcore-11.3.0                       libgeotiff/1.7.1-GCCcore-11.3.0
   ICU/71.1-GCCcore-11.3.0                             libgit2/1.4.3-GCCcore-11.3.0
   IPython/8.5.0-GCCcore-11.3.0                        libglvnd/1.4.0-GCCcore-11.3.0
   ImageMagick/7.1.0-37-GCCcore-11.3.0                 libiconv/1.17-GCCcore-11.3.0
   JasPer/2.0.33-GCCcore-11.3.0                        libjpeg-turbo/2.0.6-GCCcore-11.2.0
   Java/11.0.20                                (11)    libjpeg-turbo/2.1.3-GCCcore-11.3.0  (D)
   Judy/1.0.5-GCCcore-11.3.0                           libogg/1.3.5-GCCcore-11.3.0
   LAME/3.100-GCCcore-11.3.0                           libopus/1.3.1-GCCcore-11.3.0
   LLVM/14.0.3-GCCcore-11.3.0                          libpciaccess/0.16-GCCcore-11.3.0
   LZO/2.10-GCCcore-11.3.0                             libpng/1.6.37-GCCcore-11.2.0
   LibTIFF/4.3.0-GCCcore-11.3.0                        libpng/1.6.37-GCCcore-11.3.0        (D)
   LittleCMS/2.13.1-GCCcore-11.3.0                     libreadline/8.1-GCCcore-11.2.0
   M4/1.4.18                                           libreadline/8.1.2-GCCcore-11.3.0
   M4/1.4.19-GCCcore-11.2.0                            libreadline/8.2-GCCcore-12.2.0      (D)
   M4/1.4.19-GCCcore-11.3.0                            libsndfile/1.1.0-GCCcore-11.3.0
   M4/1.4.19-GCCcore-12.2.0                            libsodium/1.0.18-GCCcore-11.3.0
   M4/1.4.19                                   (D)     libtirpc/1.3.2-GCCcore-11.3.0
   MPFR/4.1.0-GCCcore-11.3.0                           libtool/2.4.6-GCCcore-11.2.0
   Mako/1.2.0-GCCcore-11.3.0                           libtool/2.4.7-GCCcore-11.3.0
   Maven/3.6.3                                         libtool/2.4.7                       (D)
   Mesa/22.0.3-GCCcore-11.3.0                          libunwind/1.6.2-GCCcore-11.3.0
   Meson/0.62.1-GCCcore-11.3.0                         libvorbis/1.3.7-GCCcore-11.3.0
   MetaGeneAnnotator/20080819-x86-64                   libxml2/2.9.10-GCCcore-11.2.0
   NASM/2.15.05-GCCcore-11.2.0                         libxml2/2.9.13-GCCcore-11.3.0       (D)
   NASM/2.15.05-GCCcore-11.3.0                 (D)     libxslt/1.1.34-GCCcore-11.3.0
   NLopt/2.7.1-GCCcore-11.3.0                          libyaml/0.2.5-GCCcore-11.3.0
   Ninja/1.10.2-GCCcore-11.3.0                         lxml/4.9.1-GCCcore-11.3.0
   OpenBLAS/0.3.20-GCC-11.3.0                          lz4/1.9.3-GCCcore-11.3.0
   OpenJPEG/2.5.0-GCCcore-11.3.0                       make/4.3-GCCcore-11.3.0
   OpenMPI/4.1.4-GCC-11.3.0                            ncurses/6.2-GCCcore-11.2.0
   OpenPGM/5.2.122-GCCcore-11.3.0                      ncurses/6.2
   OpenSSL/1.1                                         ncurses/6.3-GCCcore-11.3.0
   PCRE/8.45-GCCcore-11.2.0                            ncurses/6.3-GCCcore-12.2.0          (D)
   PCRE/8.45-GCCcore-11.3.0                    (D)     netCDF/4.9.0-gompi-2022a
   PCRE2/10.40-GCCcore-11.3.0                          nettle/3.8-GCCcore-11.3.0
   PMIx/4.1.2-GCCcore-11.3.0                           networkx/2.8.4-foss-2022a
   PROJ/9.0.0-GCCcore-11.3.0                           nlohmann_json/3.10.5-GCCcore-11.3.0
   Pango/1.50.7-GCCcore-11.3.0                         nodejs/16.15.1-GCCcore-11.3.0
   Perl/5.34.0-GCCcore-11.2.0                          numactl/2.0.14-GCCcore-11.3.0
   Perl/5.34.1-GCCcore-11.3.0                  (D)     pixman/0.40.0-GCCcore-11.3.0
   PyTables/3.8.0-foss-2022a                           pkg-config/0.29.2-GCCcore-11.2.0
   PyYAML/6.0-GCCcore-11.3.0                           pkg-config/0.29.2-GCCcore-11.3.0    (D)
   Python/3.9.6-GCCcore-11.2.0-bare                    pkgconf/1.8.0-GCCcore-11.3.0
   Python/3.10.4-GCCcore-11.3.0-bare                   pkgconf/1.8.0
   Python/3.10.4-GCCcore-11.3.0                (D)     pkgconf/1.9.3-GCCcore-12.2.0        (D)
   R/4.2.1-foss-2022a                                  py-cpuinfo/9.0.0-GCCcore-11.3.0
   RE2/2022-06-01-GCCcore-11.3.0                       pybind11/2.9.2-GCCcore-11.3.0
   RapidJSON/1.1.0-GCCcore-11.3.0                      ruamel.yaml/0.17.21-GCCcore-11.3.0
   Rust/1.60.0-GCCcore-11.3.0                          scikit-build/0.15.0-GCCcore-11.3.0
   SQLite/3.36-GCCcore-11.2.0                          scikit-learn/1.1.2-foss-2022a
   SQLite/3.38.3-GCCcore-11.3.0                (D)     snakemake/7.22.0-foss-2022a
   SWIG/4.0.2-GCCcore-11.2.0                           snappy/1.1.9-GCCcore-11.3.0
   ScaLAPACK/2.2.0-gompi-2022a-fb                      tabixpp/1.1.2-GCC-11.3.0
   SciPy-bundle/2022.05-foss-2022a                     tantan/40-GCC-11.2.0
   Szip/2.1.1-GCCcore-11.3.0                           tbl2asn/20230713-linux64
   Tcl/8.6.11-GCCcore-11.2.0                           unimap/0.1-GCCcore-11.3.0
   Tcl/8.6.12-GCCcore-11.3.0                   (D)     utf8proc/2.7.0-GCCcore-11.3.0
   Tk/8.6.12-GCCcore-11.3.0                            util-linux/2.37-GCCcore-11.2.0
   UCC/1.0.0-GCCcore-11.3.0                            util-linux/2.38-GCCcore-11.3.0      (D)
   UCX/1.12.1-GCCcore-11.3.0                           vispr/0.4.14-foss-2022a
   UDUNITS/2.2.28-GCCcore-11.3.0                       webin-cli/1.8.9
   USEARCH/11.0.667-i86linux32                         wgsim/20111017-GCC-11.2.0
   UnZip/6.0-GCCcore-11.2.0                            xorg-macros/1.19.3-GCCcore-11.3.0
   UnZip/6.0-GCCcore-11.3.0                    (D)     zlib/1.2.11-GCCcore-11.2.0
   UniFrac/1.3.2-foss-2022a                            zlib/1.2.11
   VSEARCH/2.22.1-GCC-11.3.0                           zlib/1.2.12-GCCcore-11.3.0
   VarScan/2.4.4-Java-11                               zlib/1.2.12-GCCcore-12.2.0
   VirSorter2/2.2.4-foss-2022a                         zlib/1.2.12                         (D)
   WFA2/2.3.3-GCCcore-11.3.0                           zstd/1.5.2-GCCcore-11.3.0
   X11/20220504-GCCcore-11.3.0

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
dp@r03:~$ cat /etc/os-release
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
[dp@node-08-1 ~]$ cat /etc/os-release
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
