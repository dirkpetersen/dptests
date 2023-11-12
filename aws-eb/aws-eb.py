#! /usr/bin/env python3

"""
AWS-EB builds Easybuild packages on AWS EC2 instances 
and uploads them to S3 buckets for later use.
"""
# internal modules
import sys, os, argparse, json, configparser, tarfile 
import urllib3, datetime, tarfile, zipfile, textwrap, platform  
import hashlib, math, signal, shlex, time, re, inspect
import shutil, tempfile, glob, subprocess, socket, traceback
if sys.platform.startswith('linux'):
    import getpass, pwd, grp
# stuff from pypi
import requests, boto3, botocore, psutil, packaging
try:
    from easybuild.framework.easyconfig.parser import EasyConfigParser
    from easybuild.tools.build_log import EasyBuildError    
except:
    print('Error: EasyBuild not found. Please install it first.')

__app__ = 'AWS-EB, a user friendly build tool for AWS EC2'
__version__ = '0.1.0.33'

def main():
        
    if args.debug:
        pass

    if len(sys.argv) == 1:        
        print(textwrap.dedent(f'''\n
            For example, use one of these commands:
              aws-eb config 
              aws-eb launch
              aws-eb download
              aws-eb ssh
            '''))

    # Instantiate classes required by all functions         
    cfg = ConfigManager(args)
    bld = Builder(args, cfg)
    aws = AWSBoto(args, cfg, bld)

    if args.version:
        args_version(cfg)


    # call a function for each sub command in our CLI
    if args.subcmd in ['config', 'cnf']:
        subcmd_config(args, cfg, aws)
    elif args.subcmd in ['launch', 'lau']:
        subcmd_launch(args, cfg, bld, aws)
    elif args.subcmd in ['download', 'bld']:
        subcmd_download(args, cfg, bld, aws)
    elif args.subcmd in ['ssh', 'scp']: #or args.unmount:
        subcmd_ssh(args, cfg, aws)

def args_version(cfg):
    print(f'AWS-EB version: {__version__}')
    print(f'Python version:\n{sys.version}')
    try:
        print('Rclone version:', subprocess.run([os.path.join(cfg.binfolderx, 'rclone'), '--version'], 
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True).stdout.split('\n')[0])
    except FileNotFoundError as e:
        print(f'Error: {e}')
        return False
    return True
    
def subcmd_config(args, cfg, aws):
    # configure user and / or team settings 
    # arguments are Class instances passed from main

    first_time=True
    
    if not cfg.binfolder:
        binfolder = '~/.local/bin'
        cfg.binfolderx = os.path.expanduser(binfolder)
        if not os.path.exists(cfg.binfolderx):
            os.makedirs(cfg.binfolderx, mode=0o775)
        cfg.write('general', 'binfolder', binfolder)
    else:
        first_time=False

    urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
    if not cfg.read('general', 'no-rclone-download'):
        
        if os.path.exists(os.path.join(cfg.binfolderx,'rclone')):
            if os.path.exists(os.path.join(cfg.binfolderx,'bak.rclone')):
                os.remove(os.path.join(cfg.binfolderx,'bak.rclone'))
            os.rename(os.path.join(cfg.binfolderx,'rclone'),os.path.join(cfg.binfolderx,'bak.rclone'))
        print(" Installing rclone ... please wait ... ", end='', flush=True)
        if platform.machine() in ['arm64', 'aarch64']:
            rclone_url = 'https://downloads.rclone.org/rclone-current-linux-arm64.zip'
        else:
            rclone_url = 'https://downloads.rclone.org/rclone-current-linux-amd64.zip'
        cfg.copy_binary_from_zip_url(rclone_url, 'rclone', 
                            '/rclone-v*/',cfg.binfolderx)
        print("Done!",flush=True)

    # general setup 
    defdom = cfg.get_domain_name()
    whoami = getpass.getuser()

    if args.monitor:
        # monitoring only setup, do not continue 
        me = os.path.join(cfg.binfolderx,'aws-eb.py')
        cfg.write('general', 'email', args.monitor)
        cfg.add_systemd_cron_job(f'{me} launch --monitor','30')
        return True
    
    print('\n*** Asking a few questions ***')
    print('*** For most you can just hit <Enter> to accept the default. ***\n')

    # set the correct permission for cfg.config_root 
    try:
        os.chmod(cfg.config_root, 0o2775)
    except:
        pass

    # domain-name not needed right now
    #domain = cfg.prompt('Enter your domain name:',
    #                    f'{defdom}|general|domain','string')
    emailaddr = cfg.prompt('Enter your email address:',
                            f'{whoami}@{defdom}|general|email','string')
    emailstr = emailaddr.replace('@','-')
    emailstr = emailstr.replace('.','-')


    print("")

    # cloud setup
    bucket = cfg.prompt('Please confirm/edit S3 bucket name to be created in all used profiles.',
                        f'aws-eb-{emailstr}|general|bucket','string')
    archiveroot = cfg.prompt('Please confirm/edit the root path inside your S3 bucket',
                                'archive|general|archiveroot','string')
    s3_storage_class =  cfg.prompt('Please confirm/edit the AWS S3 Storage class',
                                'INTELLIGENT_TIERING|general|s3_storage_class','string')


    # if there is a shared ~/.aws/config copy it over
    if cfg.config_root_local != cfg.config_root:
        cfg.replicate_ini('ALL',cfg.awsconfigfileshr,cfg.awsconfigfile)
    cfg.create_aws_configs()

    aws_region = cfg.get_aws_region('aws')
    if not aws_region:
        aws_region = cfg.get_aws_region()

    if not aws_region:
        aws_region =  cfg.prompt('Please select AWS S3 region (e.g. us-west-2 for Oregon)',
                                aws.get_aws_regions())
    aws_region =  cfg.prompt('Please confirm/edit the AWS S3 region', aws_region)
            
    #cfg.create_aws_configs(None, None, aws_region)
    print(f"\n  Verify that bucket '{bucket}' is configured ... ")
    
    allowed_aws_profiles = ['default', 'aws', 'AWS'] # for accessing glacier use one of these    
    profs = cfg.get_aws_profiles()

    for prof in profs:
        if prof in allowed_aws_profiles:
            cfg.set_aws_config(prof, 'region', aws_region)
            if prof == 'AWS' or prof == 'aws':
                cfg.write('general', 'aws_profile', prof)                    
            elif prof == 'default': 
                cfg.write('general', 'aws_profile', 'default')
            aws.create_s3_bucket(bucket, prof)
        
    print('\nDone!\n')


def subcmd_launch(args,cfg,bld,aws):

    cfg.printdbg ("build:", args.awsprofile)
    cfg.printdbg(f'default cmdline: aws-eb build')

    if args.monitor:
        # aws inactivity and cost monitoring
        aws.monitor_ec2()
        return True

    if args.awsprofile and args.awsprofile not in cfg.get_aws_profiles():
        print(f'Profile "{args.awsprofile}" not found.')
        return False    
    
    if args.list:
        # list all folders in the archive
        print('\nAll available EC2 instance families:')
        print("--------------------------------------")
        fams = aws.get_ec2_instance_families()
        print(' '.join(fams))        
        print("\nGPU    Instance Families")
        print("--------------------------")
        for c, i in aws.gpu_types.items():
            print(f'{c}: {i}')               
        print("\nCPU Type    Instance Families")
        print("--------------------------------")
        for c, i in aws.cpu_types.items():
            print(f'{c}: {" ".join(i)}')
        return True
    
    # GPU types trump CPU types 
    fams_c = []
    fam = ''
    if args.cputype:
        fams_c = aws.get_ec2_instance_families_from_cputype(args.cputype)
        if not fams_c:
            print(f'CPU type "{args.cputype}" not found. Run build --list to see types.')
            return False
        fam = fams_c[0]

    if args.gputype:
        fam = aws.get_ec2_instance_families_from_gputype(args.gputype)
        if not fam:
            print(f'GPU type "{args.gputype}" not found. Run build --list to see types.')
            return False
        args.cputype =  aws.get_ec2_cputype_from_instance_family(fam)

    if not args.cputype:
        print('Please specify a CPU or a GPU type. Run build --list to see types.')
        return False

    os_id, version_id = cfg.get_os_release_info()

    if not os_id or not version_id:
        print('Could not determine OS release information.')
        return False    
    
    s3_prefix = f'{os_id}-{version_id}_{args.cputype}'
    if args.gputype:
        s3_prefix += f'_{args.gputype}'

    instance_type = aws.get_ec2_cheapest_instance_type(fam, args.vcpus, args.mem*1024)
    print('Cheapest:', instance_type)

    if not args.build:
        aws.ec2_deploy(256, instance_type) # 256GB disk for the build instance
        return True

    # *******************************************
    # Start EasyBuild process here 
    print('s3_prefix:', s3_prefix)
    ecfgroot = os.path.join(cfg.home_dir, '.local', 'easybuild', 'easyconfigs')
    bld.build_all(ecfgroot, s3_prefix, bio_only=args.bioonly)

    #if not aws.check_bucket_access_folders(args.folders):
    #    return False
        
def subcmd_download(args,cfg,bld,aws):

    cfg.printdbg ("restore:",args.cores, args.awsprofile, args.noslurm, 
        args.days, args.retrieveopt, args.nodownload, args.folders)
    fld = '" "'.join(args.folders)
    cfg.printdbg(f'default cmdline: aws-eb restore "{fld}"')

    # ********* 
    if args.monitor:
        # aws inactivity and cost monitoring
        aws.monitor_ec2()
        return True
    
    if args.awsprofile and args.awsprofile not in cfg.get_aws_profiles():
        print(f'Profile "{args.awsprofile}" not found.')
        return False    
    if not aws.check_bucket_access_folders(args.folders):
        return False
    
    if args.ec2:
        # run ec2_deploy(self, bucket='', prefix='', recursive=False, profile=None):
        ret = aws.ec2_deploy(args.folders)
        return True


def subcmd_ssh(args, cfg, aws):

    ilist = aws.ec2_list_instances('Name', 'AWSEBSelfDestruct')
    ips = [sublist[0] for sublist in ilist if sublist]
    if args.list:
        if ips:
            print("Running EC2 Instances:")
            for row in ilist:
                print(' - '.join(row))        
        else:
            print('No running instances detected')
        return True        
    if args.terminate:
        aws.ec2_terminate_instance(args.terminate)
        return True        
    myhost = cfg.read('cloud', 'ec2_last_instance')        
    if ips and not myhost in ips:
        print(f'{myhost} is no longer running, replacing with {ips[-1]}')
        myhost = ips[-1]
        #cfg.write('cloud', 'ec2_last_instance', myhost)
    if args.subcmd == 'ssh':
        if args.sshargs:
            myhost = args.sshargs[0]
        print(f'Connecting to {myhost} ...')
        aws.ssh_execute('ec2-user', myhost)
        return True
    elif args.subcmd == 'scp':
        if len(args.sshargs) != 2:
            print('The "scp" sub command supports currently 2 arguments')
            return False
        hostloc = next((i for i, item in enumerate(args.sshargs) if ":" in item), None)
        if hostloc == 0:
            # the hostname is in the first argument: download
            host, remote_path = args.sshargs[0].split(':')
            ret=aws.ssh_download('ec2-user', host, remote_path, args.sshargs[1])
        elif hostloc == 1:
            # the hostname is in the second argument: uploaad
            host, remote_path = args.sshargs[1].split(':')
            ret=aws.ssh_upload('ec2-user', host, args.sshargs[0], remote_path)
        else:
            print('The "scp" sub command supports currently 2 arguments')
            return False
        print(ret.stdout,ret.stderr)

class Builder:
    def __init__(self, args, cfg):        
        self.args = args
        self.cfg = cfg
        self.rclone_download_compare = '--checksum'
        self.rclone_upload_compare = '--checksum'
        self.min_toolchains = {'system': 'system', 'GCC': '11.0', 'GCCcore' : '11.0', 'LLVM' : '12.0', 'foss' : '2023a'}
        #self.min_toolchains = {'system': 'system', 'GCC': '11.0', 'GCCcore': '11.0', 'foss': '2022a', 'fosscuda': '2021a'}
        self.eb_root = os.path.join('/', 'opt', 'eb')

    def build_all(self, easyconfigroot, s3_prefix, bio_only=False):

        # install a lot of required junk 
        self._install_os_dependencies(easyconfigroot)
        softwaredir = os.path.join(self.eb_root, 'software')

        # build all new easyconfigs in a folder tree
        for root, dirs, files in self._walker(easyconfigroot):
            print(f'  Processing folder "{root}" newest easyconfigs... ')
            try:
                if easyconfigroot==root:
                    # main directory do something there
                    pass
                ebfile = self._get_latest_easyconfig(root)
                if not ebfile:
                    continue
                print(f'  Processing easyconfig "{ebfile}" ... ', flush=True)
                ebpath = os.path.join(root, ebfile)
                if not os.path.isfile(ebpath):
                    continue 
                name, version, tc, dep, cls, instdir = self._read_easyconfig(ebpath)
                if name in self.min_toolchains.keys(): # if this is the toolchain package itself    
                    if self.cfg.sversion(version) < self.cfg.sversion(self.min_toolchains[name]):
                        print(f'  * Easyconfig {name} version {version} too old.', flush=True)
                        continue
                if tc['name'] not in self.min_toolchains.keys():
                    print(f'  * Toolchain {tc["name"]} not supported.', flush=True)
                    continue
                if self.cfg.sversion(tc['version']) < self.cfg.sversion(self.min_toolchains[tc['name']]):
                    print(f'  * Toolchain version {tc["version"]} of {tc["name"]} too old.', flush=True)
                    continue
                if cls != 'bio' and bio_only:
                    # we want to may be only build bio packages
                    print(f'  {name} is not a bio package', flush=True)
                    continue
                if dep:
                    print(f'  installing OS dependencies: {dep}', flush=True)
                    self._install_packages(dep)
                # install easybuild package 
                print(f" Downloading previous packages ... ", flush=True)
                self.download(f':s3:{self.cfg.archivepath}', self.eb_root, s3_prefix)
                print(f" Unpacking previous packages ... ", flush=True)
                all_tars, new_tars = self._untar_eb_software(softwaredir)
                print(f" Installing {ebfile} ... ", flush=True)
                ret = subprocess.run(['eb', '--robot', '--umask=002', ebpath])
                print(f'*** EASYBUILD RETURNCODE: {ret.returncode}', flush=True)
                print(f" Tarring up new packages ... ", flush=True)
                all_tars, new_tars = self._tar_eb_software(softwaredir)
                print(f" Uploading new packages ... ", flush=True)
                self.upload(self.eb_root, f':s3:{self.cfg.archivepath}', s3_prefix)
                                                
            except subprocess.CalledProcessError:                
                print(f"  Builder.build_all: A CalledProcessError occurred while building {ebfile}.", flush=True)
                ## make sure we store the logfile
                continue

            except Exception as e:
                print(f"  Builder.build_all: An unexpected error occurred:\n{e}", flush=True)
                traceback.print_exc()
                continue
        try:
            print(f" Final upload using checksums ... ", flush=True)
            self.rclone_upload_compare = '--checksum'
            self.upload(self.eb_root, f':s3:{self.cfg.archivepath}', s3_prefix)
        except Exception as e:
            print(f"  Builder.build_all(final): An unexpected error occurred when uploading:\n{e}", flush=True)
            pass

        return True
    
    def _install_os_dependencies(self, easyconfigroot):
        # install OS dependencies from all easyconfigs (~ 400 packages)        
        package_skip_set = set() # avoid duplicates
        self._install_packages(['pigz'], package_skip_set)        
        for root, dirs, files in self._walker(easyconfigroot):
            print(f'  Processing folder "{root}" for OS depts... ')
            for ebfile in files:
                if ebfile.endswith('.eb'):
                    ebpath = os.path.join(root, ebfile)
                    _, _, _, dep, _, _ = self._read_easyconfig(ebpath)
                    if dep:
                        print(f'  installing OS dependencies: {dep}')
                        self._install_packages(dep, package_skip_set)
                        for package_tuple in dep: # avoid duplicates
                            if isinstance(package_tuple, str):
                                package_tuple = (package_tuple,)                            
                            for package_name in package_tuple:
                                package_skip_set.add(package_name)

    def _tar_folder_old(self, folder):
        # Ensure the directory exists
        if not os.path.isdir(folder):
            raise ValueError(f"The directory {folder} does not exist.")    
        # Define the name of the tarball
        tarball_name = f'{folder}.tar.gz'
        # Create a tar.gz archive
        with tarfile.open(tarball_name, 'w:gz') as tar:
            # Add the directory to the tarball
            tar.add(folder, arcname=os.path.basename(folder))
        print(f'Directory {folder} has been archived as {tarball_name}')

    def _tar_eb_software(self, folder):
        new_tars = []
        all_tars = []
        for root, dirs, files in self._walker(folder):
            # Check if 'easybuild' is in the directories
            if 'easybuild' in dirs:
                # Extract the folder name which should be the version, and the parent folder which should be the package
                version_dir = os.path.basename(root)
                package_dir = os.path.basename(os.path.dirname(root))
                package_root = os.path.dirname(root)

                # Create the tarball name
                tarball_name = f'{package_dir}-{version_dir}.eb.tar.gz'
                tarball_path = os.path.join(folder, package_dir, tarball_name)
                all_tars.append(tarball_path)
                if os.path.isfile(tarball_path):
                    print(f"Tarball {tarball_path} already exists ...")
                    continue
                new_tars.append(tarball_path)

                # Print info for the user
                print(f"Creating tarball {tarball_path} from {root}...", flush=True)

                # Use tar with pigz for compression
                try:
                    subprocess.run([
                        "tar",
                        "-I", f"pigz -p {self.args.vcpus}",  # Call pigz for compression with X CPUs
                        "-cf", f'{tarball_path}.tmp',  # Create and verbosely list files processed
                        "-C", package_root,  # Change to the parent directory of version
                        version_dir  # Specify the directory to compress
                    ], check=True)
                    os.rename(f'{tarball_path}.tmp', tarball_path)
                    print(f"Successfully created tarball: {tarball_path}")
                except subprocess.CalledProcessError as e:
                    print(f"An error occurred while creating tarball: {e}")
        return all_tars, new_tars
    
    def _untar_eb_software(self, folder):
        new_tars = []
        all_tars = []
        for root, dirs, files in self._walker(folder):
            # Extract package name from the root directory
            package_name = os.path.basename(root)
            for filename in files:
                if filename.endswith('.eb.tar.gz'):
                    # Strip the '.tar.gz' extension and then extract the version
                    version = filename.replace('.eb.tar.gz', '').replace(package_name + '-', '')

                    # Construct the expected path for the version directory
                    version_dir_path = os.path.join(root, version)

                    # Check if the 'easybuild' directory exists within the version directory
                    easybuild_path = os.path.join(version_dir_path, 'easybuild')
                    file_path = os.path.join(root, filename)
                    all_tars.append(file_path)
                    if not os.path.exists(easybuild_path):
                        print(f"Unpacking {file_path} into {version_dir_path}...", flush=True)
                        # try:
                        #     # Decompress with pigz through tar command
                        #     subprocess.run([
                        #         "tar",
                        #         "-I", f"pigz -p {self.args.vcpus}",
                        #         "-xf", file_path,
                        #         "-C", root 
                        #     ], check=True)
                        #     print(f"Successfully unpacked: {file_path}")
                        #     new_tars.append(file_path)
                        # except subprocess.CalledProcessError as e:
                        #     print(f"An error occurred while unpacking {file_path}: {e}")
                        try:
                            # Extract the tar.gz file using tarfile
                            with tarfile.open(file_path, "r:gz") as tar:
                                tar.extractall(path=root)
                            print(f"Successfully unpacked: {file_path}")
                        except Exception as e:
                            print(f"An error occurred while unpacking {file_path}: {e}")

                    else:
                        pass
                        #print(f"Skipping unpacking of {file_path} as 'easybuild' directory already exists in {version_dir_path}.")
        return all_tars, new_tars

    def _get_latest_easyconfig(self,directory):
 
        from packaging.version import parse, InvalidVersion
        version_file_dict = {}
        version_pattern = re.compile(r'-(\d+(?:\.\d+)*)(?:-(\w+(?:-\d+(?:\.\d+)*(?:[ab]\d+)?)?))?\.')

        for filename in os.listdir(directory):
            match = version_pattern.search(filename)
            if match:
                software_version = match.group(1)
                toolchain_version = match.group(2) if match.group(2) else '0'  # Default to '0' if no toolchain
                # Attempt to parse the toolchain version as a semantic version
                try:
                    toolchain_version_parsed = parse(toolchain_version)
                except InvalidVersion:
                    # If it's not a valid semantic version, we will use the string itself for comparison
                    toolchain_version_parsed = toolchain_version
                
                version_tuple = (parse(software_version), toolchain_version_parsed)
                version_file_dict[version_tuple] = filename
        
        if not version_file_dict:
            return None

        # Sort by software version and then toolchain version
        # Custom sort to handle non-standard version strings
        def sort_key(version_tuple):
            software, toolchain = version_tuple
            if isinstance(toolchain, str):
                # Assume non-standard versions are older, set them as the minimum
                return software, parse('0')
            return software, toolchain

        latest_version = sorted(version_file_dict.keys(), key=sort_key, reverse=True)[0]
        return version_file_dict[latest_version]

    def _read_easyconfig(self, ebfile):
        
        ec_dict = {}
        try:
            # Initialize EasyConfigParser with the easyconfig file
            ec_dict = EasyConfigParser(ebfile).get_config_dict()
        except EasyBuildError as e:
            print("An error occurred while parsing the easyconfig file:", e)

        toolchain = ec_dict.get('toolchain', {})
        name = ec_dict.get('name', "")
        version = ec_dict.get('version', "")
        versionsuffix = ec_dict.get('versionsuffix', "")

        toolchain_str = f"-{toolchain['name']}-{toolchain['version']}" if toolchain['name'] != 'system' else ""

        # Construct the version suffix string
        version_suffix_str = f"{version}{versionsuffix}" if versionsuffix else version

        # Construct the installation directory path
        install_dir = f"{name}/{version_suffix_str}{toolchain_str}"

        return name, version, toolchain, ec_dict.get('osdependencies', ""), ec_dict.get('moduleclass', ""), install_dir

    def _get_os_type(self):
        os_info = {}
        if os.path.isfile("/etc/os-release"):
            with open("/etc/os-release") as f:
                for line in f:
                    key, value = line.strip().split("=", 1)
                    os_info[key] = value.strip('"')
        # Prioritize ID_LIKE over ID
        os_type = os_info.get('ID_LIKE', os_info.get('ID', None))
        # Handle the case where ID_LIKE can contain multiple space-separated strings
        if os_type and ' ' in os_type:
            os_type = os_type.split(' ')[0]  # Get the first 'like' identifier
        return os_type

    def _install_packages(self, os_dependencies, package_skip_set=[]):
        os_type = self._get_os_type()        
        # Determine the appropriate package manager for the detected OS type
        package_manager = None
        if os_type in ['debian', 'ubuntu']:
            package_manager = 'apt'
        elif os_type in ['fedora', 'centos', 'redhat', 'rhel']:
            package_manager = 'dnf'
        
        if not package_manager:
            print("Unsupported operating system.")
            return
        
        for package_tuple in os_dependencies:
            if isinstance(package_tuple, str):
                package_tuple = (package_tuple,)
            installed = False
            for package_name in package_tuple:
                # Check if the package has a known OS-specific suffix
                if package_name in package_skip_set:
                    print(f"Skipping {package_name} because it was already installed.")
                    continue
                if (package_name.endswith('-dev') and os_type in ['debian', 'ubuntu']) or \
                (package_name.endswith('-devel') and os_type in ['fedora', 'centos', 'redhat']):
                    try:
                        print(f"Installing {package_name} with {package_manager}")                    
                        subprocess.run(['sudo', package_manager, 'install', '-y', package_name], check=True)
                        installed = True
                        break
                    except subprocess.CalledProcessError:
                        pass
            
            # If none of the packages in the tuple had a suffix, we try to install each until one succeeds
            if not installed:
                for package_name in package_tuple:
                    if package_name in package_skip_set:
                        print(f"Skipping {package_name} because it was already installed.")
                        continue                    
                    try:
                        print(f"Attempting to install {package_name} with {package_manager}")
                        subprocess.run(['sudo', package_manager, 'install', '-y', package_name], check=True)
                        print(f"Installed {package_name} successfully.")
                        break  # Stop trying after the first successful install
                    except subprocess.CalledProcessError:
                        # If the package installation failed, it might be the wrong package for the OS,
                        # so continue trying the next packages in the tuple
                        pass

    def upload(self, source, target, s3_prefix):

        source = os.path.abspath(source)
    
        rclone = Rclone(self.args,self.cfg)

        if not self.rclone_upload_compare == '--size-only':
            print ('  Uploading Bootstrap output ... ', flush=True)
            ret = rclone.copy(os.path.expanduser('~/'),
                            f'{target}/{s3_prefix}/logs/',
                            '--include', 'out.bootstrap.*' 
                            )        

        print ('  Uploading Modules ... ', flush=True)
        ret = rclone.copy(os.path.join(source,'modules'),
                          f'{target}/{s3_prefix}/modules/', 
                          '--links', self.rclone_upload_compare
                        )

        print ('  Uploading Sources ... ', flush=True)
        ret = rclone.copy(os.path.join(source,'sources'),
                          f'{target}/sources/', 
                          '--links', self.rclone_upload_compare
                        )

        print ('  Uploading Software ... ', flush=True)
        ret = rclone.copy(os.path.join(source,'software'),
                          f'{target}/{s3_prefix}/software/', 
                          '--links', self.rclone_upload_compare,
                          '--include', '*.eb.tar.gz' 
                        )
        
        print ('  Uploading EB output ... ', flush=True)
        ret = rclone.copy(os.path.expanduser('~/'),
                          f'{target}/{s3_prefix}/logs/',
                          '--include', 'out.easybuild.*' 
                        )        

        # after the first successful upload do a size only compare
        self.rclone_upload_compare  = '--size-only'
                
        self.cfg.printdbg('*** RCLONE copy ret ***:\n', ret, '\n')
        #print ('Message:', ret['msg'].replace('\n',';'))
        if ret['stats']['errors'] > 0:
            print('Last Error:', ret['stats']['lastError'])
            print('Copying was not successful.')
            return False
                
        ttransfers=ret['stats']['totalTransfers']
        tbytes=ret['stats']['totalBytes']
        total=self.convert_size(tbytes)
        if self.args.debug:
            print('\n')
            print('Speed:', ret['stats']['speed'])
            print('Transfers:', ret['stats']['transfers'])
            print('Tot Transfers:', ret['stats']['totalTransfers'])
            print('Tot Bytes:', ret['stats']['totalBytes'])
            print('Tot Checks:', ret['stats']['totalChecks'])

        #   {'bytes': 0, 'checks': 0, 'deletedDirs': 0, 'deletes': 0, 'elapsedTime': 2.783003019, 
        #    'errors': 1, 'eta': None, 'fatalError': False, 'lastError': 'directory not found', 
        #    'renames': 0, 'retryError': True, 'speed': 0, 'totalBytes': 0, 'totalChecks': 0, 
        #    'totalTransfers': 0, 'transferTime': 0, 'transfers': 0}   
        
        print(f'Upload finished. {ttransfers} files with {total} transferred.\n')

    def download(self, source, target, s3_prefix=None):
               
        rclone = Rclone(self.args,self.cfg)
            
        print ('  Downloading Modules ... ', flush=True)
        ret = rclone.copy(f'{source}/{s3_prefix}/modules/',
                          os.path.join(target,'modules'), 
                          '--links', self.rclone_download_compare
                        )

        print ('  Downloading Sources ... ', flush=True)
        ret = rclone.copy(f'{source}/sources/',
                          os.path.join(target,'sources'), 
                          '--links', self.rclone_download_compare
                        )
        
        self._make_files_executable(os.path.join(target,'sources','generic'))

        print ('  Downloading Software ... ', flush=True)
        ret = rclone.copy(f'{source}/{s3_prefix}/software/',
                          os.path.join(target,'software'), 
                          '--links', self.rclone_download_compare, 
                          '--include', '*.eb.tar.gz' 
                        )
        
        # for subsequent download comparison size is enough
        self.rclone_download_compare = '--size-only'

        self.cfg.printdbg('*** RCLONE copy ret ***:\n', ret, '\n')
        #print ('Message:', ret['msg'].replace('\n',';'))
        if ret['stats']['errors'] > 0:
            print('Last Error:', ret['stats']['lastError'])
            print('Copying was not successful.')
            return False
            # lastError could contain: Object in GLACIER, restore first

        ttransfers=ret['stats']['totalTransfers']
        tbytes=ret['stats']['totalBytes']
        total=self.convert_size(tbytes)
        if self.args.debug:
            print('\n')
            print('Speed:', ret['stats']['speed'])
            print('Transfers:', ret['stats']['transfers'])
            print('Tot Transfers:', ret['stats']['totalTransfers'])
            print('Tot Bytes:', ret['stats']['totalBytes'])
            print('Tot Checks:', ret['stats']['totalChecks'])

        #   {'bytes': 0, 'checks': 0, 'deletedDirs': 0, 'deletes': 0, 'elapsedTime': 2.783003019, 
        #    'errors': 1, 'eta': None, 'fatalError': False, 'lastError': 'directory not found', 
        #    'renames': 0, 'retryError': True, 'speed': 0, 'totalBytes': 0, 'totalChecks': 0, 
        #    'totalTransfers': 0, 'transferTime': 0, 'transfers': 0}   
        # checksum

 
        print(f'Download finished. {ttransfers} files with {total} transferred.')
   
        return -1

    def _make_files_executable(self, path):
        for root, dirs, files in self._walker(path):
            for file in files:
                if not file.endswith('.tar.gz'):
                    file_path = os.path.join(root, file)
                    if not os.access(file_path, os.X_OK):
                        print(f'Making {file_path} executable')
                        os.chmod(file_path, os.stat(file_path).st_mode | 0o111)

    def test_write(self, directory):
        testpath=os.path.join(directory,'.aws-eb.test')
        try:
            with open(testpath, "w") as f:
                f.write('just a test')
            os.remove(testpath)
            return True
        except PermissionError as e:
            if e.errno == 13:  # Check if error number is 13 (Permission denied)
                #print("Permission denied. Please ensure you have the necessary permissions to access the file or directory.")
                return 13
            else:
                print(f"An unexpected PermissionError occurred in {directory}:\n{e}")            
                return False
        except Exception as e:
            if e.errno == 2:
                #No such file or directory:
                return 2
            else:
                print(f"An unexpected error occurred in {directory}:\n{e}")
                return False

        

    def _get_file_stats(self,filepath):
        try:
            # Use lstat to get stats of symlink itself, not the file it points to
            stats = os.lstat(filepath)
            return stats.st_size, stats.st_mtime, stats.st_atime
        except FileNotFoundError:
            print(f"{filepath} not found.")
            return None, None, None
            
    
    
    def md5sumex(self, file_path):
        try:
            cmd = f'md5sum {file_path}'
            ret = subprocess.run(cmd, stdout=subprocess.PIPE, 
                            stderr=subprocess.PIPE, Shell=True)                    
            if ret.returncode != 0:
                print(f'md5sum return code > 0: {cmd} Error:\n{ret.stderr}')
            return ret.stdout.strip() #, ret.stderr.strip()

        except Exception as e:
            print (f'md5sum Error: {str(e)}')
            return None, str(e)
                             
    def md5sum(self, file_path):
        md5_hash = hashlib.md5()
        with open(file_path, "rb") as f:
            for chunk in iter(lambda: f.read(4096), b""):
                md5_hash.update(chunk)
        return md5_hash.hexdigest()

    
    def uid2user(self,uid):
        # try to convert uid to user name
        try:
            return pwd.getpwuid(uid)[0]
        except:
            self.cfg.printdbg(f'uid2user: Error converting uid {uid}')
            return uid

    def gid2group(self,gid):
        # try to convert gid to group name
        try:
            return grp.getgrgid(gid)[0]
        except:
            self.cfg.printdbg(f'gid2group: Error converting gid {gid}')
            return gid

    def daysago(self,unixtime):
        # how many days ago is this epoch time ?
        if not unixtime: 
            self.cfg.printdbg('daysago: an integer is required (got type NoneType)')
            return 0
        diff=datetime.datetime.now()-datetime.datetime.fromtimestamp(unixtime)
        return diff.days
    
    def convert_size(self, size_bytes):
        if size_bytes == 0:
            return "0B"
        size_name = ("B", "KiB", "MiB", "GiB", "TiB", "PiB", "EiB")
        i = int(math.floor(math.log(size_bytes, 1024)))
        p = math.pow(1024, i)
        s = round(size_bytes/p, 3)
        return f"{s} {size_name[i]}"    
    
    def _get_newest_file_atime(self, folder_path, folder_atime=None):
        # Because the folder atime is reset when crawling we need
        # to lookup the atime of the last accessed file in this folder
        if not os.path.exists(folder_path) or not os.path.isdir(folder_path):
            if self.args.debug and not self.args.pwalkcsv:
                print(f" Invalid folder path: {folder_path}")
            return folder_atime
        last_accessed_time = None
        #last_accessed_file = None
        for file_name in os.listdir(folder_path):
            file_path = os.path.join(folder_path, file_name)
            if os.path.isfile(file_path):
                accessed_time = os.path.getatime(file_path)
                if last_accessed_time is None or accessed_time > last_accessed_time:
                    last_accessed_time = accessed_time
                    #last_accessed_file = file_path
        if last_accessed_time == None:
            last_accessed_time = folder_atime
        return last_accessed_time
    

    def _walker(self, top, skipdirs=['.snapshot', '__archive__']):
        """ returns subset of os.walk  """
        for root, dirs, files in os.walk(top,topdown=True,onerror=self._walkerr): 
            for skipdir in skipdirs:
                if skipdir in dirs:
                    dirs.remove(skipdir)  # don't visit this directory 
            yield root, dirs, files 

    def _walkerr(self, oserr):    
        sys.stderr.write(str(oserr))
        sys.stderr.write('\n')
        return 0

    def _get_last_directory(self, path):
        # Remove any trailing slashes
        path = path.rstrip(os.path.sep)
        # Split the path by the separator
        path_parts = path.split(os.path.sep)
        # Return the last directory
        return path_parts[-1]
 
    def _get_mount_info(self,fs_types=None):
        file_path='/proc/self/mountinfo'
        if fs_types is None:
            fs_types = {'nfs', 'nfs4', 'cifs', 'smb', 'afs', 'ncp', 
                        'ncpfs', 'glusterfs', 'ceph', 'beegfs', 
                        'lustre', 'orangefs', 'wekafs', 'gpfs'}
        mountinfo_list = []
        with open(file_path, 'r') as f:
            for line in f:
                fields = line.strip().split(' ')
                _, _, _, _, mount_point, _ = fields[:6]
                for field in fields[6:]:
                    if field == '-':
                        break
                fs_type, mount_source, _ = fields[-3:]
                mount_source_folder = mount_source.split(':')[-1] if ':' in mount_source else ''
                if fs_type in fs_types:
                    mountinfo_list.append({
                        'mount_source_folder': mount_source_folder,
                        'mount_point': mount_point,
                        'fs_type': fs_type,
                        'mount_source': mount_source,
                    })
        return mountinfo_list
        

    def download_restored_file(self, bucket_name, object_key, local_path):
        s3 = boto3.resource('s3')
        s3.Bucket(bucket_name).download_file(object_key, local_path)
        print(f'Downloaded {object_key} to {local_path}.')

    def _upload_file_to_s3(self, filename, bucket, object_name=None, profile=None):
        session = boto3.Session(profile_name=profile) if profile else boto3.Session()
        s3 = session.client('s3')
        # If S3 object_name was not specified, use the filename
        if object_name is None:
            object_name = os.path.basename(filename)
        try:
            # Upload the file with Intelligent-Tiering storage class
            s3.upload_file(filename, bucket, object_name, ExtraArgs={'StorageClass': 'INTELLIGENT_TIERING'})
            self.cfg.printdbg(f"File {object_name} uploaded to Intelligent-Tiering storage class!")
            #print(f"File {filename} uploaded successfully to Intelligent-Tiering storage class!")
        except Exception as e:
            print(f"An error occurred: {e}")
            return False
        return True

    # initiate_restore(bucket_name, object_key, restore_days)    
    # while not check_restore_status(bucket_name, object_key):
    #     print('Waiting for restoration to complete...')
    #     time.sleep(60)  # Wait 60 seconds before checking again
    # download_restored_file(bucket_name, object_key, local_path)

class Rclone:
    def __init__(self, args, cfg):
        self.args = args
        self.cfg = cfg
        self.rc = os.path.join(self.cfg.binfolderx,'rclone')

    # ensure that file exists or nagging /home/dp/.config/rclone/rclone.conf

    #backup: rclone --verbose --files-from tmpfile --use-json-log copy --max-depth 1 ./tests/ :s3:posix-dp/tests4/ --exclude .aws-eb.md5sum
    #restore: rclone --verbose --use-json-log copy --max-depth 1 :s3:posix-dp/tests4/ ./tests2
    #rclone copy --verbose --use-json-log --max-depth 1  :s3:posix-dp/tests5/ ./tests5
    #rclone --use-json-log checksum md5 ./tests/.aws-eb.md5sum :s3:posix-dp/tests2/
    # storage tier for each file 
    #rclone lsf --csv :s3:posix-dp/tests4/ --format=pT
    # list without subdir 
    #rclone lsjson --metadata --no-mimetype --no-modtime --hash :s3:posix-dp/tests4
    #rclone checksum md5 ./tests/.aws-eb.md5sum --verbose --use-json-log :s3:posix-dp/archive/home/dp/gh/aws-eb/tests

    def _run_rc(self, command):

        command = self._add_opt(command, '--verbose')
        command = self._add_opt(command, '--use-json-log')
        self.cfg.printdbg('Rclone command:', " ".join(command))
        try:
            ret = subprocess.run(command, capture_output=True, text=True, env=self.cfg.envrn)
            if ret.returncode != 0:
                #pass
                sys.stderr.write(f'*** Error, Rclone return code > 0:\n{ret.stderr} Command:\n{" ".join(command)}\n\n')
                # list of exit codes 
                # 0 - success
                # 1 - Syntax or usage error
                # 2 - Error not otherwise categorised
                # 3 - Directory not found
                # 4 - File not found
                # 5 - Temporary error (one that more retries might fix) (Retry errors)
                # 6 - Less serious errors (like 461 errors from dropbox) (NoRetry errors)
                # 7 - Fatal error (one that more retries won't fix, like account suspended) (Fatal errors)
                # 8 - Transfer exceeded - limit set by --max-transfer reached
                # 9 - Operation successful, but no files transferred
            
            #lines = ret.stderr.decode('utf-8').splitlines() #needed if you do not use ,text=True
            #locked_dirs = '\n'.join([l for l in lines if "Locked Dir:" in l]) 
            #print("   STDOUT:",ret.stdout)
            #print("   STDERR:",ret.stderr)
            #rclone mount --daemon
            return ret.stdout.strip(), ret.stderr.strip()

        except Exception as e:
            print (f'Rclone Error: {str(e)}')
            return None, str(e)

    def _run_bk(self, command):
        #command = self._add_opt(command, '--verbose')
        #command = self._add_opt(command, '--use-json-log')
        cmdline=" ".join(command)
        self.cfg.printdbg('Rclone command:', cmdline)
        try:
            ret = subprocess.Popen(command, preexec_fn=os.setsid, stdin=subprocess.PIPE, 
                        stdout=subprocess.PIPE, text=True, env=self.cfg.envrn)
            #_, stderr = ret.communicate(timeout=3)  # This does not work with rclone
            if ret.stderr:
                sys.stderr.write(f'*** Error in command "{cmdline}":\n {ret.stderr} ')
            return ret.pid
        except Exception as e:
            print (f'Rclone Error: {str(e)}')
            return None

    def copy(self, src, dst, *args):
        command = [self.rc, 'copy'] + list(args)
        command.append(src)  #command.append(f'{src}/')
        command.append(dst)
        out, err = self._run_rc(command)
        if out:
            print(f'rclone copy output: {out}')
        #print('ret', err)
        stats, ops = self._parse_log(err) 
        if stats:
            return stats[-1] # return the stats
        else:
            return []
    
        #b'{"level":"warning","msg":"Time may be set wrong - time from \\"posix-dp.s3.us-west-2.amazonaws.com\\" is -9m17.550965814s different from this computer","source":"fshttp/http.go:200","time":"2023-04-16T14:40:47.44907-07:00"}'    

    def checksum(self, md5file, dst, *args):
        #checksum md5 ./tests/.aws-eb.md5sum
        command = [self.rc, 'checksum'] + list(args)
        command.append('md5')
        command.append(md5file)
        command.append(dst)
        #print("Command:", command)
        out, err = self._run_rc(command)
        if out:
            print(f'rclone checksum output: {out}')
        #print('ret', err)
        stats, ops = self._parse_log(err) 
        if stats:
            return stats[-1] # return the stats
        else:
            return []

    def mount(self, url, mountpoint, *args):
        if not shutil.which('fusermount3'):
            print('Could not find "fusermount3". Please install the "fuse3" OS package')
            return False
        if not url.endswith('/'): url+'/'
        mountpoint = mountpoint.rstrip(os.path.sep)
        command = [self.rc, 'mount'] + list(args)
        # might use older rclone, if fuse3 is not installed
        #if os.path.isfile('/usr/bin/rclone'):
        #    command = ['/usr/bin/rclone', 'mount'] + list(args)            
        #command.append('--daemon') # not reliable, just starting background process
        try:
            #os.chmod(mountpoint, 0o2775)
            current_permissions = os.stat(mountpoint).st_mode
            new_permissions = (current_permissions & ~0o07) | 0o05
            os.chmod(mountpoint, new_permissions)            
        except:
            pass
        command.append('--allow-non-empty')
        command.append('--default-permissions')
        command.append('--read-only')
        command.append('--no-checksum')
        command.append('--quiet')
        command.append(url)
        command.append(mountpoint)
        pid = self._run_bk(command)
        return pid

    def unmount(self, mountpoint, wait=False):
        mountpoint = mountpoint.rstrip(os.path.sep)
        if self._is_mounted(mountpoint):
            if shutil.which('fusermount3'):
                cmd = ['fusermount3', '-u', mountpoint]
                ret = subprocess.run(cmd, capture_output=False, text=True, env=self.cfg.envrn)
            else:
                rclone_pids = self._get_pids('rclone')
                fld_pids = self._get_pids(mountpoint, True)
                common_pids = [value for value in rclone_pids if value in fld_pids]
                for pid in common_pids:
                    try:
                        os.kill(pid, signal.SIGTERM)
                        if wait:
                            _, _ = os.waitpid(int(pid), 0)
                        return True
                    except PermissionError:
                        print(f'Permission denied when trying to send signal SIGTERM to rclone process with PID {pid}.')
                    except Exception as e:
                        print(f'An unexpected error occurred when trying to send signal SIGTERM to rclone process with PID {pid}: {e}') 
        else:
            print(f'\nError: Folder {mountpoint} is currently not used as a mountpoint by rclone.')
                
    def version(self):
        command = [self.rc, 'version']
        return self._run_rc(command)

    def get_mounts(self):
        mounts = []
        with open('/proc/mounts', 'r') as f:
            for line in f:
                parts = line.split()
                mount_point, fs_type = parts[1], parts[2]
                if fs_type.startswith('fuse.rclone'):
                    mounts.append(mount_point)
        return mounts

    def _get_pids(self, process, full=False):
        process = process.rstrip(os.path.sep)
        if full:
            command = ['pgrep', '-f', process]
        else:
            command = ['pgrep', process]
        try:
            output = subprocess.check_output(command)
            pids = [int(pid) for pid in output.decode().split('\n') if pid]
            return pids
        except subprocess.CalledProcessError:
            # No rclone processes found
            return []

    def _is_mounted(self, folder_path):
        folder_path = os.path.realpath(folder_path)  # Resolve any symbolic links
        with open('/proc/mounts', 'r') as f:
            for line in f:
                parts = line.split()
                mount_point, fs_type = parts[1], parts[2]
                if mount_point == folder_path and fs_type.startswith('fuse.rclone'):
                    return True


    def _add_opt(self, cmd, option, value=None):
        if option in cmd:
            return cmd
        cmd.append(option)
        if value:
            cmd.append(value)
        return cmd
    
    def _parse_log(self, strstderr):
        lines=strstderr.split('\n')
        data = [json.loads(line.rstrip()) for line in lines if line[0] == "{"]
        stats = []
        operations = []
        for obj in data:
            if 'accounting/stats' in obj['source']:
                stats.append(obj)
            elif 'operations/operations' in obj['source']:
                operations.append(obj)
        return stats, operations

        # stats":{"bytes":0,"checks":0,"deletedDirs":0,"deletes":0,"elapsedTime":4.121489785,"errors":12,"eta":null,"fatalError":false,
        # "lastError":"failed to open source object: Object in GLACIER, restore first: bucket=\"posix-dp\", key=\"tests4/table_example.py\"",
        # "renames":0,"retryError":true,"speed":0,"totalBytes":0,"totalChecks":0,"totalTransfers":0,"transferTime":0,"transfers":0},
        # "time":"2023-04-16T10:18:46.121921-07:00"}


class AWSBoto:
    # we write all config entries as files to '~/.config'
    # to make it easier for bash users to read entries 
    # with a simple var=$(cat ~/.config/aws-eb/section/entry)
    # entries can be strings, lists that are written as 
    # multi-line files and dictionaries which are written to json

    def __init__(self, args, cfg, bld):
        self.args = args
        self.cfg = cfg
        self.bld = bld
        self.awsprofile = self.cfg.awsprofile
        self.cpu_types = {
            "graviton-2": ['m6g','c6g', 'c6gn', 't4g' ,'g5g'],
            "graviton-3": ['m7g', 'c7g', 'c7gn'],
            "epyc-gen-1": ['t3a'],
            "epyc-gen-2": ['m5a', 'c5a', 'r5a', 'g4ad', 'p4', 'inf2', 'g5'],
            "epyc-gen-3": ['m6a', 'c6a', 'r6a', 'p5'],
            "epyc-gen-4": ['m7a', 'c7a', 'r7a'],
            "xeon-gen-1": ['m4', 'c4', 't2', 'r4', 'p3' ,'p2', 'f1', 'g3'],
            "xeon-gen-2": ['m5', 'c5', 'c5n', 'm5n', 'm5zn', 'r5', 't3', 't3n', 'dl1', 'inf1', 'g4dn', 'vt1'],
            "xeon-gen-3": ['m6i', 'c6i', 'm6in', 'c6in', 'r6i', 'r6id', 'r6idn', 'r6in', 'trn1'],
            "xeon-gen-4": ['m7i', 'c7i', 'm7i-flex'],
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


    def get_ec2_instance_families_from_cputype(self, cpu_type):
        return self.cpu_types.get(cpu_type,[])

    def get_ec2_instance_families_from_gputype(self, gpu_type):
        return self.gpu_types.get(gpu_type,"")
    
    def get_ec2_cputype_from_instance_family(self, ifamily):
        for cputype, families in self.cpu_types.items():
            if ifamily in families:
                return cputype
        return ""

    def get_ec2_instance_families(self, profile=None):
        session = boto3.Session(profile_name=profile) if profile else boto3.Session()
        ec2 = session.client('ec2')
        families = set()
        try:
            paginator = ec2.get_paginator('describe_instance_types')
            for page in paginator.paginate():
                for itype in page['InstanceTypes']:
                    # Extract the family (prefix before the dot) and add it to the set
                    family = itype['InstanceType'].split('.')[0]
                    families.add(family)

        except Exception as e:
            print(f"Error retrieving instance types: {e}")
            return

        # Convert the set to a list and sort it to list families in order
        sorted_families = sorted(list(families))
        return sorted_families
    

    def get_ec2_cheapest_instance_type(self, family, min_vcpu, min_memory, gpu_type=None, profile=None):
        session = boto3.Session(profile_name=profile) if profile else boto3.Session()
        ec2 = session.client('ec2')

        # Initialize variables
        suitable_types = []
        try:
            paginator = ec2.get_paginator('describe_instance_types')
            for page in paginator.paginate():
                for itype in page['InstanceTypes']:
                    # Check if the instance type belongs to the specified family
                    if itype['InstanceType'].startswith(family):
                        vcpus = itype['VCpuInfo']['DefaultVCpus']
                        memory = itype['MemoryInfo']['SizeInMiB']

                        # Check if the instance meets the minimum vCPU and memory requirements
                        if vcpus >= min_vcpu and memory >= min_memory:
                            suitable_types.append(itype)

            # Sort suitable types by vCPUs and memory to try to get the smallest (and possibly cheapest) type
            suitable_types.sort(key=lambda x: (x['VCpuInfo']['DefaultVCpus'], x['MemoryInfo']['SizeInMiB']))

            # Assuming the first instance type is the cheapest based on the sorting
            if suitable_types:
                return suitable_types[0]['InstanceType']
            else:
                return "No suitable instance type found."

        except Exception as e:
            print(f"Error retrieving instance types: {e}")
            return None

    def get_aws_regions(self, profile=None, provider='AWS'):
        # returns a list of AWS regions 
        if provider == 'AWS':
            try:
                session = boto3.Session(profile_name=profile) if profile else boto3.Session()
                regions = session.get_available_regions('ec2')
                # make the list a little shorter 
                regions = [i for i in regions if not i.startswith('ap-')]
                return sorted(regions, reverse=True)
            except:
                return ['us-west-2','us-west-1', 'us-east-1', '']

    # def check_bucket_access_folders(self, folders, readwrite=False):
    #     # check all the buckets that have been used for archiving 
    #     sufficient = True
    #     myaccess = 'read'
    #     if readwrite:
    #         myaccess = 'write'
    #     buckets = []
    #     for folder in folders:
    #         bucket, *_ = self.bld.archive_get_bucket_info(folder)
    #         buckets.append(bucket)
    #     buckets = list(set(buckets)) # remove dups
    #     for bucket in buckets:
    #         if not self.check_bucket_access(bucket, readwrite):
    #             print (f' You have no {myaccess} access to bucket "{bucket}" !')
    #             sufficient = False
    #     return sufficient

    def check_bucket_access(self, bucket_name, readwrite=False, profile=None):
        
        if not bucket_name:
            print('check_bucket_access: bucket_name empty. You may have not yet configured a S3 bucket name. Please run "aws-eb config" first')
            sys.exit(1)    
        if not self._check_s3_credentials(profile):
            print('_check_s3_credentials failed. Please edit file ~/.aws/credentials')
            return False
        session = boto3.Session(profile_name=profile) if profile else boto3.Session()
        ep_url = self.cfg._get_aws_s3_session_endpoint_url(profile)
        s3 = session.client('s3', endpoint_url=ep_url)
        
        try:
            # Check if bucket exists
            s3.head_bucket(Bucket=bucket_name)
        except botocore.exceptions.ClientError as e:
            error_code = e.response['Error']['Code']
            if error_code == '403':
                print(f"Error: Access denied to bucket {bucket_name} for profile {self.awsprofile}. Check your permissions.")
            elif error_code == '404':
                print(f"Error: Bucket {bucket_name} does not exist in profile {self.awsprofile}.")
                print("run 'aws-eb config' to create this bucket.")
            else:
                print(f"Error accessing bucket {bucket_name} in profile {self.awsprofile}: {e}")
            return False
        except Exception as e:
            print(f"An unexpected error in function check_bucket_access for profile {self.awsprofile}: {e}")
            return False

        if not readwrite:
            return True
        
        # Test write access by uploading a small test file
        try:
            test_object_key = "test_write_access.txt"
            s3.put_object(Bucket=bucket_name, Key=test_object_key, Body="Test write access")
            #print(f"Successfully wrote test to {bucket_name}")

            # Clean up by deleting the test object
            s3.delete_object(Bucket=bucket_name, Key=test_object_key)
            #print(f"Successfully deleted test object from {bucket_name}")
            return True
        except botocore.exceptions.ClientError as e:
            print(f"Error: cannot write to bucket {bucket_name} in profile {self.awsprofile}: {e}")
            return False
        
        
####


    def create_s3_bucket(self, bucket_name, profile=None):   
        if not self._check_s3_credentials(profile, verbose=True):
            print(f"Cannot create bucket '{bucket_name}' with these credentials")
            print('check_s3_credentials failed. Please edit file ~/.aws/credentials')
            return False 
        region = self.cfg.get_aws_region(profile)
        session = boto3.Session(profile_name=profile) if profile else boto3.Session()
        ep_url = self.cfg._get_aws_s3_session_endpoint_url(profile)
        s3_client = session.client('s3', endpoint_url=ep_url)        
        existing_buckets = s3_client.list_buckets()
        for bucket in existing_buckets['Buckets']:
            if bucket['Name'] == bucket_name:
                self.cfg.printdbg(f'S3 bucket {bucket_name} exists')
                return True
        try:
            if region and region != 'default-placement':
                response = s3_client.create_bucket(
                    Bucket=bucket_name,
                    CreateBucketConfiguration={'LocationConstraint': region}
                    )
            else:
                response = s3_client.create_bucket(
                    Bucket=bucket_name,
                    )              
            print(f"Created S3 Bucket '{bucket_name}'")
        except botocore.exceptions.BotoCoreError as e:
            print(f"BotoCoreError: {e}")
        except botocore.exceptions.ClientError as e:
            error_code = e.response['Error']['Code']
            if error_code == 'InvalidBucketName':
                print(f"Error: Invalid bucket name '{bucket_name}'\n{e}")
            elif error_code == 'BucketAlreadyExists':
                pass
                #print(f"Error: Bucket '{bucket_name}' already exists.")
            elif error_code == 'BucketAlreadyOwnedByYou':
                pass
                #print(f"Error: You already own a bucket named '{bucket_name}'.")
            elif error_code == 'InvalidAccessKeyId':
                #pass
                print("Error: InvalidAccessKeyId. The AWS Access Key Id you provided does not exist in our records")
            elif error_code == 'SignatureDoesNotMatch':
                pass
                #print("Error: Invalid AWS Secret Access Key.")
            elif error_code == 'AccessDenied':
                print("Error: Access denied. Check your account permissions for creating S3 buckets")
            elif error_code == 'IllegalLocationConstraintException':
                print(f"Error: The specified region '{region}' is not valid.")
            else:
                print(f"ClientError: {e}")
            return False
        except Exception as e:            
            print(f"An unexpected error occurred: {e}")
            return False
        encryption_configuration = {
            'Rules': [
                {
                    'ApplyServerSideEncryptionByDefault': {
                        'SSEAlgorithm': 'AES256'
                    }
                }
            ]
        }
        try:
            response = s3_client.put_bucket_encryption(
                Bucket=bucket_name,
                ServerSideEncryptionConfiguration=encryption_configuration
            )            
            print(f"Applied AES256 encryption to S3 bucket '{bucket_name}'")    
        except botocore.exceptions.ClientError as e:
            error_code = e.response['Error']['Code']
            if error_code == 'InvalidBucketName':
                print(f"Error: Invalid bucket name '{bucket_name}'\n{e}")
            elif error_code == 'AccessDenied':
                print("Error: Access denied. Check your account permissions for creating S3 buckets")
            elif error_code == 'IllegalLocationConstraintException':
                print(f"Error: The specified region '{region}' is not valid.")
            elif error_code == 'InvalidLocationConstraint':
                if not ep_url:
                    # do not show this error with non AWS endpoints 
                    print(f"Error: The specified location-constraint '{region}' is not valid")
            else:
                print(f"ClientError: {e}")                        
        except Exception as e:            
            print(f"An unexpected error occurred in create_s3_bucket: {e}")
            return False            
        return True
    
    def _check_s3_credentials(self, profile=None, verbose=False):
        session = boto3.Session(profile_name=profile) if profile else boto3.Session()
        try:
            if verbose or self.args.debug:
                self.cfg.printdbg(f'  Checking credentials for profile "{profile}" ... ', end='')            
            ep_url = self.cfg._get_aws_s3_session_endpoint_url(profile)
            s3_client = session.client('s3', endpoint_url=ep_url)            
            s3_client.list_buckets()
            if verbose or self.args.debug:
                pass
                #print('Done.')                
        #return True
        #try:
        #    pass
        except botocore.exceptions.NoCredentialsError:
            print("No AWS credentials found. Please check your access key and secret key.")
        except botocore.exceptions.EndpointConnectionError:
            print("Unable to connect to the AWS S3 endpoint. Please check your internet connection.")
        except botocore.exceptions.ClientError as e:
            error_code = e.response.get('Error', {}).get('Code')
            #error_code = e.response['Error']['Code']             
            if error_code == 'RequestTimeTooSkewed':
                print(f"The time difference between S3 storage and your computer is too high:\n{e}")
            elif error_code == 'InvalidAccessKeyId':                
                print(f"Error: Invalid AWS Access Key ID in profile {profile}:\n{e}")
                print(f"Fix your credentials in ~/.aws/credentials for profile {profile}")
            elif error_code == 'SignatureDoesNotMatch':                
                if "Signature expired" in str(e): 
                    print(f"Error: Signature expired. The system time of your computer is likely wrong:\n{e}")
                    return False
                else:
                    print(f"Error: Invalid AWS Secret Access Key in profile {profile}:\n{e}")         
            elif error_code == 'InvalidClientTokenId':
                print(f"Error: Invalid AWS Access Key ID or Secret Access Key !")
                print(f"Fix your credentials in ~/.aws/credentials for profile {profile}")                
            else:
                print(f"Error validating credentials for profile {profile}: {e}")
                print(f"Fix your credentials in ~/.aws/credentials for profile {profile}")
            return False
        except Exception as e:
            print(f"An unexpected Error in _check_s3_credentials with profile {profile}: {e}")
            sys.exit(1)
        return True
    
    # def _get_s3_data_size(self, folders, profile=None):
    #     """
    #     Get the size of data in GiB aggregated from multiple 
    #     S3 buckets from aws-eb archives identified by a 
    #     list of folders 

    #     :return: Size of the data in GiB.
    #     """
    #     session = boto3.Session(profile_name=profile) if profile else boto3.Session()
    #     s3 = session.client('s3')

    #     # Initialize total size
    #     total_size_bytes = 0
        
    #     #bucket_name, prefix, recursive=False
    #     for fld in folders:
    #         buc, pre, recur, _ = self.arch.archive_get_bucket_info(fld)
    #         # returns bucket(str), prefix(str), recursive(bool), glacier(bool)
    #         # Use paginator to handle buckets with large number of objects
    #         paginator = s3.get_paginator('list_objects_v2')
    #         for page in paginator.paginate(Bucket=buc, Prefix=pre):
    #             if "Contents" in page:  # Ensure there are objects under the specified prefix
    #                 for obj in page['Contents']:
    #                     key = obj['Key']
    #                     if recur or (key.count('/') == pre.count('/') and key.startswith(pre)):
    #                         total_size_bytes += obj['Size']
        
    #     total_size_gib = total_size_bytes / (1024 ** 3)  # Convert bytes to GiB
    #     return total_size_gib

    def ec2_deploy(self, disk_gib, instance_type, awsprofile=None):

        if not awsprofile: 
            awsprofile = self.cfg.awsprofile
        prof = self._ec2_create_iam_policy_roles_ec2profile()            
        iid, ip = self._ec2_create_instance(disk_gib, instance_type, prof, awsprofile)
        print(' Waiting for ssh host to become ready ...')
        if not self.cfg.wait_for_ssh_ready(ip):
            return False

        bootstrap_build = self._ec2_user_space_script(iid)        

        ### this block may need to be moved to a function
        argl = ['--ec2', '-e']
        cmdlist = [item for item in sys.argv if item not in argl]
        argl = ['--instance-type', '-i'] # if found remove option and next arg
        cmdlist = [x for i, x in enumerate(cmdlist) if x \
                   not in argl and (i == 0 or cmdlist[i-1] not in argl)]
        if not '--profile' in cmdlist and self.args.awsprofile:
            cmdlist.insert(1,'--profile')
            cmdlist.insert(2, self.args.awsprofile)
        if not '--build' in cmdlist:
            cmdlist.append('--build')
        cmdline = 'aws-eb.py ' + " ".join(map(shlex.quote, cmdlist[1:])) #original cmdline
        ### end block 

        print(f" will execute '{cmdline}' on {ip} ... ")
        bootstrap_build += '\n' + cmdline + f' > ~/out.easybuild.{ip}.txt 2>&1'        
        # once everything is done, commit suicide, but only if ~/no-terminate does not exist:
        bootstrap_build += f'\n[ ! -f ~/no-terminate ] && aws-eb.py ssh --terminate {iid}'
        ret = self.ssh_upload('ec2-user', ip,
            self._ec2_easybuildrc(), "easybuildrc", is_string=True)
        ret = self.ssh_upload('ec2-user', ip,
            bootstrap_build, "bootstrap.sh", is_string=True)        
        if ret.stdout or ret.stderr:
            print(ret.stdout, ret.stderr)
        ret = self.ssh_execute('ec2-user', ip, 
            'mkdir -p ~/.config/aws-eb/general')
        if ret.stdout or ret.stderr:
            print(ret.stdout, ret.stderr)        
        ret = self.ssh_upload('ec2-user', ip,
            "~/.config/aws-eb/general/*", ".config/aws-eb/general/")
        if ret.stdout or ret.stderr:
            print(ret.stdout, ret.stderr)        
        ret = self.ssh_execute('ec2-user', ip, 
            f'nohup bash bootstrap.sh < /dev/null > out.bootstrap.{ip}.txt 2>&1 &')
        if ret.stdout or ret.stderr:
            print(ret.stdout, ret.stderr)
        print(' Executed bootstrap and build script ... you may have to wait a while ...')
        print(' but you can already login using "aws-eb ssh"')

        os.system(f'echo "touch ~/no-terminate" >> ~/.bash_history.tmp')
        os.system(f'echo "grep -A1 ^ERROR: ~/out.easybuild.{ip}.txt" >> ~/.bash_history.tmp')
        os.system(f'echo "tail -f ~/out.easybuild.{ip}.txt" >> ~/.bash_history.tmp')
        os.system(f'echo "tail -f ~/out.bootstrap.{ip}.txt" >> ~/.bash_history.tmp')
        ret = self.ssh_upload('ec2-user', ip,
            "~/.bash_history.tmp", ".bash_history")
        if ret.stdout or ret.stderr:
            print(ret.stdout, ret.stderr)

        self.send_email_ses('', '', 'AWS-EB restore on EC2', f'this command line was executed on host {ip}:\n{cmdline}')

    def _ec2_create_or_get_iam_policy(self, pol_name, pol_doc, profile=None):
        session = boto3.Session(profile_name=profile) if profile else boto3.Session()
        iam = session.client('iam')

        policy_arn = None
        try:
            response = iam.create_policy(
                PolicyName=pol_name,
                PolicyDocument=json.dumps(pol_doc)
            )
            policy_arn = response['Policy']['Arn']
            print(f"Policy created with ARN: {policy_arn}")
        except iam.exceptions.EntityAlreadyExistsException as e:
            policies = iam.list_policies(Scope='Local')  
               # Scope='Local' for customer-managed policies, 
               # 'AWS' for AWS-managed policies            
            for policy in policies['Policies']:
                if policy['PolicyName'] == pol_name:
                    policy_arn = policy['Arn']
                    break
            print(f'Policy {pol_name} already exists')
        except botocore.exceptions.ClientError as e:
            error_code = e.response['Error']['Code']
            if error_code == 'AccessDenied':
                self.cfg.printdbg(f'Access denied! Please check your IAM permissions. \n   Error: {e}')
            else:
                print(f'Client Error: {e}')
        except Exception as e:
            print('Other Error:', e)
        return policy_arn

    def _ec2_create_aws_eb_iam_policy(self, profile=None):
        # Initialize session with specified profile or default
        session = boto3.Session(profile_name=profile) if profile else boto3.Session()

        # Create IAM client
        iam = session.client('iam')

        # Define policy name and policy document
        policy_name = 'AWS-EBEC2DescribePolicy'
        policy_document = {
            "Version": "2012-10-17",
            "Statement": [
                {
                    "Effect": "Allow",
                    "Action": "ec2:Describe*",
                    "Resource": "*"
                }
            ]
        }

        # Get current IAM user's details
        user = iam.get_user()
        user_name = user['User']['UserName']

        # Check if policy already exists for the user
        existing_policies = iam.list_user_policies(UserName=user_name)
        if policy_name in existing_policies['PolicyNames']:
            print(f"{policy_name} already exists for user {user_name}.")
            return

        # Create policy for user
        iam.put_user_policy(
            UserName=user_name,
            PolicyName=policy_name,
            PolicyDocument=json.dumps(policy_document)
        )

        print(f"Policy {policy_name} attached successfully to user {user_name}.")


    def _ec2_create_iam_policy_roles_ec2profile(self, profile=None):
        # create all the IAM requirement to allow an ec2 instance to
        # 1. self destruct, 2. monitor cost with CE and 3. send emails via SES
        session = boto3.Session(profile_name=profile) if profile else boto3.Session()
        iam = session.client('iam')

      # Step 0: Create IAM self destruct and EC2 read policy 
        policy_document = {
            "Version": "2012-10-17",
            "Statement": [
                {
                    "Effect": "Allow",
                    "Action": [
                        "ec2:Describe*",     # Basic EC2 read permissions
                    ],
                    "Resource": "*"
                },                
                {
                    "Effect": "Allow",
                    "Action": "ec2:TerminateInstances",
                    "Resource": "*",
                    "Condition": {
                        "StringEquals": {
                            "ec2:ResourceTag/Name": "AWSEBSelfDestruct"
                        }
                    }
                },
                {
                    "Effect": "Allow",
                    "Action": [
                        "ce:GetCostAndUsage"
                    ],
                    "Resource": "*"
                }
            ]
        }
        policy_name = 'AWSEBSelfDestructPolicy'     

        destruct_policy_arn = self._ec2_create_or_get_iam_policy(
            policy_name, policy_document, profile)
    
        # 1. Create an IAM role
        trust_policy = {
            "Version": "2012-10-17",
            "Statement": [
                {
                    "Effect": "Allow",
                    "Principal": {"Service": "ec2.amazonaws.com"},
                    "Action": "sts:AssumeRole"
                },
            ]
        }

        role_name = "AWS-EBEC2Role"
        try:
            iam.create_role(
                RoleName=role_name,
                AssumeRolePolicyDocument=json.dumps(trust_policy),
                Description='AWS-EB role allows Billing, SES and Terminate'
            )
        except iam.exceptions.EntityAlreadyExistsException:        
            print (f'Role {role_name} already exists.') 
        except botocore.exceptions.ClientError as e:
            error_code = e.response['Error']['Code']
            if error_code == 'AccessDenied':
                self.cfg.printdbg(f'Access denied! Please check your IAM permissions. \n   Error: {e}')
            else:
                print(f'Client Error: {e}')
        except Exception as e:            
            print('Other Error:', e)
        
        # 2. Attach permissions policies to the IAM role
        cost_explorer_policy = "arn:aws:iam::aws:policy/AWSBillingReadOnlyAccess"
        ses_policy = "arn:aws:iam::aws:policy/AmazonSESFullAccess"

        try:
        
            iam.attach_role_policy(
                RoleName=role_name,
                PolicyArn=cost_explorer_policy
            )
            
            iam.attach_role_policy(
                RoleName=role_name,
                PolicyArn=ses_policy
            )

            iam.attach_role_policy(
                RoleName=role_name,
                PolicyArn=destruct_policy_arn
            )
        except iam.exceptions.PolicyNotAttachableException as e:
            print(f"Policy {e.policy_arn} is not attachable. Please check your permissions.")
            return False
        except botocore.exceptions.ClientError as e:
            error_code = e.response['Error']['Code']
            if error_code == 'AccessDenied':
                self.cfg.printdbg(f'Access denied! Please check your IAM permissions. \n   Error: {e}')
            else:
                print(f'Client Error: {e}')
        except Exception as e:
            print('Other Error:', e)
            return False
        # 3. Create an instance profile and associate it with the role
        instance_profile_name = "AWS-EBEC2Profile"
        try:
            iam.create_instance_profile(
                InstanceProfileName=instance_profile_name
            )
            iam.add_role_to_instance_profile(
                InstanceProfileName=instance_profile_name,
                RoleName=role_name
            )
        except iam.exceptions.EntityAlreadyExistsException:
            print (f'Profile {instance_profile_name} already exists.')
            return instance_profile_name
        except botocore.exceptions.ClientError as e:
            error_code = e.response['Error']['Code']
            if error_code == 'AccessDenied':
                self.cfg.printdbg(f'Access denied! Please check your IAM permissions. \n   Error: {e}')
            else:
                print(f'Client Error: {e}')
            return None
        except Exception as e:            
            print('Other Error:', e)
            return None
        
        # Give AWS a moment to propagate the changes
        print('wait for 15 sec ...')
        time.sleep(15)  # Wait for 15 seconds

        return instance_profile_name
    
    def _ec2_create_and_attach_security_group(self, instance_id, profile=None):
        session = boto3.Session(profile_name=profile) if profile else boto3.Session()
        ec2 = session.resource('ec2')
        client = session.client('ec2')

        group_name = 'SSH-HTTP-ICMP'
        
        # Check if security group already exists
        security_groups = client.describe_security_groups(Filters=[{'Name': 'group-name', 'Values': [group_name]}])
        if security_groups['SecurityGroups']:
            security_group_id = security_groups['SecurityGroups'][0]['GroupId']
        else:
            # Create security group
            response = client.create_security_group(
                GroupName=group_name,
                Description='Allows SSH and ICMP inbound traffic'
            )
            security_group_id = response['GroupId']
        
            # Allow ports 22, 80, 443, 8000-9000, ICMP
            client.authorize_security_group_ingress(
                GroupId=security_group_id,
                IpPermissions=[
                    {
                        'IpProtocol': 'tcp',
                        'FromPort': 22,
                        'ToPort': 22,
                        'IpRanges': [{'CidrIp': '0.0.0.0/0'}]
                    },
                    {
                        'IpProtocol': 'tcp',
                        'FromPort': 80,
                        'ToPort': 80,
                        'IpRanges': [{'CidrIp': '0.0.0.0/0'}]
                    },
                    {
                        'IpProtocol': 'tcp',
                        'FromPort': 443,
                        'ToPort': 443,
                        'IpRanges': [{'CidrIp': '0.0.0.0/0'}]
                    },
                    {
                        'IpProtocol': 'tcp',
                        'FromPort': 8000,
                        'ToPort': 9000,
                        'IpRanges': [{'CidrIp': '0.0.0.0/0'}]
                    },                    {
                        'IpProtocol': 'icmp',
                        'FromPort': -1,  # -1 allows all ICMP types
                        'ToPort': -1,
                        'IpRanges': [{'CidrIp': '0.0.0.0/0'}]
                    }
                ]
            )
        
        # Attach the security group to the instance
        instance = ec2.Instance(instance_id)
        current_security_groups = [sg['GroupId'] for sg in instance.security_groups]
        
        # Check if the security group is already attached to the instance
        if security_group_id not in current_security_groups:
            current_security_groups.append(security_group_id)
            instance.modify_attribute(Groups=current_security_groups)

        return security_group_id

    def _ec2_get_latest_amazon_linux2_ami(self, profile=None):
        session = boto3.Session(profile_name=profile) if profile else boto3.Session()
        ec2_client = session.client('ec2')

        response = ec2_client.describe_images(
            Filters=[
                {'Name': 'name', 'Values': ['al2023-ami-*']},
                {'Name': 'state', 'Values': ['available']},
                {'Name': 'architecture', 'Values': ['x86_64']},
                {'Name': 'virtualization-type', 'Values': ['hvm']}
            ],
            Owners=['amazon']

            #amzn2-ami-hvm-2.0.*-x86_64-gp2
            #al2023-ami-kernel-default-x86_64

        )

        # Sort images by creation date to get the latest
        images = sorted(response['Images'], key=lambda k: k['CreationDate'], reverse=True)
        if images:
            return images[0]['ImageId']
        else:
            return None        

    def _create_progress_bar(self, max_value):
        def show_progress_bar(iteration):
            percent = ("{0:.1f}").format(100 * (iteration / float(max_value)))
            length = 50  # adjust as needed for the bar length
            filled_length = int(length * iteration // max_value)
            bar = "█" * filled_length + '-' * (length - filled_length)
            print(f'\r|{bar}| {percent}%', end='\r')
            if iteration == max_value: 
                print()

        return show_progress_bar

    def _ec2_cloud_init_script(self):
        # Define the User Data script
        long_timezone = self.cfg.get_time_zone()
        userdata = textwrap.dedent(f'''
        #! /bin/bash
        dnf install -y gcc mdadm
        #bigdisks=$(lsblk --fs --json | jq -r '.blockdevices[] | select(.children == null and .fstype == null) | "/dev/" + .name')
        bigdisks='/dev/sdm'
        numdisk=$(echo $bigdisks | wc -w)
        mkdir /opt/eb
        if [[ $numdisk -gt 1 ]]; then
          mdadm --create /dev/md0 --level=0 --raid-devices=$numdisk $bigdisks
          mkfs -t xfs /dev/md0
          mount /dev/md0 /opt/eb
        elif [[ $numdisk -eq 1 ]]; then
          mkfs -t xfs $bigdisks
          mount $bigdisks /opt/eb
        fi
        chown ec2-user /opt/eb
        dnf check-update
        dnf update -y                                   
        dnf install -y at gcc vim wget python3-pip python3-psutil 
        hostnamectl set-hostname aws-eb
        timedatectl set-timezone '{long_timezone}'
        loginctl enable-linger ec2-user
        systemctl start atd
        dnf upgrade
        dnf install -y mc git docker lua lua-posix lua-devel tcl-devel nodejs-npm
        dnf group install -y 'Development Tools'
        cd /tmp
        wget https://sourceforge.net/projects/lmod/files/Lmod-8.7.tar.bz2
        tar -xjf Lmod-8.7.tar.bz2
        cd Lmod-8.7 && ./configure && make install
        ''').strip()
        return userdata
    
    def _ec2_easybuildrc(self, bscript='~/easybuildrc'):
        threads = self.args.vcpus*2
        return textwrap.dedent(f'''        
        test -d /usr/local/lmod/lmod/init && source /usr/local/lmod/lmod/init/bash
        export MODULEPATH=/opt/eb/modules/all:/opt/eb/modules/lang:/opt/eb/modules/compiler:/opt/eb/modules/ai:/opt/eb/modules/bio
        # export MODULEPATH=/opt/eb/modules/tools:/opt/eb/modules/lang:/opt/eb/modules/compiler:/opt/eb/modules/bio
        #
        export EASYBUILD_JOB_CORES={self.args.vcpus}
        export EASYBUILD_CUDA_COMPUTE_CAPABILITIES=7.5,8.0,8.6,9.0
        # export EASYBUILD_BUILDPATH=/dev/shm/$USER # could run out of space
        export EASYBUILD_PREFIX=/opt/eb
        export EASYBUILD_JOB_OUTPUT_DIR=$EASYBUILD_PREFIX/batch-output
        export EASYBUILD_DEPRECATED=5.0
        export EASYBUILD_JOB_BACKEND=Slurm
        export EASYBUILD_PARALLEL={threads}
        # export EASYBUILD_GITHUB_USER=$USER
        export  EASYBUILD_UPDATE_MODULES_TOOL_CACHE=True
        export  EASYBUILD_ROBOT_PATHS=/home/ec2-user/.local/easybuild/easyconfigs:/opt/eb/fh/fh_easyconfigs
        ''').strip()
            
    def _ec2_user_space_script(self, instance_id='', bscript='~/bootstrap.sh'):
        # Define script that will be installed by ec2-user 
        emailaddr = self.cfg.read('general','email')
        #short_timezone = datetime.datetime.now().astimezone().tzinfo
        long_timezone = self.cfg.get_time_zone()
        return textwrap.dedent(f'''
        #! /bin/bash
        mkdir -p ~/.config/aws-eb
        echo 'PS1="\\u@aws-eb:\\w$ "' >> ~/.bashrc
        echo 'source ~/easybuildrc' >> ~/.bashrc
        echo '#export EC2_INSTANCE_ID={instance_id}' >> ~/.bashrc
        echo '#export AWS_DEFAULT_REGION={self.cfg.aws_region}' >> ~/.bashrc
        echo '#export TZ={long_timezone}' >> ~/.bashrc
        echo '#alias singularity="apptainer"' >> ~/.bashrc        
        aws configure set aws_access_key_id {os.environ['AWS_ACCESS_KEY_ID']}
        aws configure set aws_secret_access_key {os.environ['AWS_SECRET_ACCESS_KEY']}
        aws configure set region {self.cfg.aws_region}
        aws configure --profile {self.cfg.awsprofile} set aws_access_key_id {os.environ['AWS_ACCESS_KEY_ID']}
        aws configure --profile {self.cfg.awsprofile} set aws_secret_access_key {os.environ['AWS_SECRET_ACCESS_KEY']}
        aws configure --profile {self.cfg.awsprofile} set region {self.cfg.aws_region}
        sed -i 's/aws_access_key_id [^ ]*/aws_access_key_id /' {bscript}
        sed -i 's/aws_secret_access_key [^ ]*/aws_secret_access_key /' {bscript}
        curl -s https://raw.githubusercontent.com/apptainer/apptainer/main/tools/install-unprivileged.sh | bash -s - ~/.local
        echo '#! /bin/bash' > ~/.local/bin/get-public-ip
        echo 'ETOKEN=$(curl -sX PUT "http://169.254.169.254/latest/api/token" -H "X-aws-ec2-metadata-token-ttl-seconds: 21600")' >> ~/.local/bin/get-public-ip
        cp -f ~/.local/bin/get-public-ip ~/.local/bin/get-local-ip
        echo 'curl -sH "X-aws-ec2-metadata-token: $ETOKEN" http://169.254.169.254/latest/meta-data/public-ipv4' >> ~/.local/bin/get-public-ip
        echo 'curl -sH "X-aws-ec2-metadata-token: $ETOKEN" http://169.254.169.254/latest/meta-data/local-ipv4' >> ~/.local/bin/get-local-ip
        chmod +x ~/.local/bin/get-public-ip
        chmod +x ~/.local/bin/get-local-ip
        # curl https://raw.githubusercontent.com/dirkpetersen/aws-eb/main/install.sh | bash 
        curl -Ls https://raw.githubusercontent.com/dirkpetersen/dptests/main/aws-eb/aws-eb.py -o ~/.local/bin/aws-eb.py
        chmod +x ~/.local/bin/aws-eb.py
        # wait for pip3 and lmod to be installed
        until [ -f /usr/bin/pip3 ]; do sleep 5; done; echo "pip3 exists, please wait ..."  
        until [ -f /usr/local/lmod/lmod/init/bash ]; do sleep 5; done; echo "lmod exists, please wait ..."  
        python3 -m pip install --upgrade wheel
        python3 -m pip install boto3 easybuild packaging
        source ~/easybuildrc
        aws-eb.py config --monitor '{emailaddr}'
        echo "" >> ~/.bash_profile
        ''').strip()
    
    def _ec2_create_instance(self, disk_gib, instance_type, iamprofile=None, profile=None):
        
        session = boto3.Session(profile_name=profile) if profile else boto3.Session()
        ec2 = session.resource('ec2')
        client = session.client('ec2')
        
        # Define the block device mapping for an EBS volume to be attached to the instance
        block_device_mappings = [
        {
            'DeviceName': '/dev/sdm',  # Ensure that this device name is supported and free in your EC2 instance
            'Ebs': {
                'VolumeSize': disk_gib,  # Volume size in GiB (1 TB = 1024 GiB)
                'DeleteOnTermination': True,  # True if the volume should be deleted after instance is terminated
                'VolumeType': 'gp2',  # The type of volume to create (gp3 is generally a good default)
                # Additional parameters like Iops and Throughput can be specified for 'gp3' volume type
            },
        }]    
             
        if self.args.instancetype:
            instance_type = self.args.instancetype 
        
        if not instance_type:
            print("No suitable instance type found!")
            return False

        # Create a new EC2 key pair
        key_path = os.path.join(self.cfg.config_root,'cloud',f'{self.cfg.ssh_key_name}.pem')
        if not os.path.exists(key_path):
            try:
                client.describe_key_pairs(KeyNames=[self.cfg.ssh_key_name])
                # If the key pair exists, delete it
                client.delete_key_pair(KeyName=self.cfg.ssh_key_name)
            except client.exceptions.ClientError:
                # Key pair doesn't exist in AWS, no need to delete
                pass                        
            key_pair = ec2.create_key_pair(KeyName=self.cfg.ssh_key_name)
            os.makedirs(os.path.join(self.cfg.config_root,'cloud'),exist_ok=True)            
            with open(key_path, 'w') as key_file:
                key_file.write(key_pair.key_material)
            os.chmod(key_path, 0o600)  # Set file permission to 600

        imageid = self._ec2_get_latest_amazon_linux2_ami(profile)
        print(f'Using Image ID: {imageid}')

        #print(f'*** userdata-script:\n{self._ec2_user_data_script()}')

        iam_instance_profile={}
        if iamprofile:
            iam_instance_profile={
                'Name': iamprofile  # Use the instance profile name
            }        
        print(f'IAM Instance profile: {iamprofile}.')    

        # iam_instance_profile = {}
        
        try:
            # Create EC2 instance
            instance = ec2.create_instances(
                ImageId=imageid,
                MinCount=1,
                MaxCount=1,
                InstanceType=instance_type,
                KeyName=self.cfg.ssh_key_name,
                UserData=self._ec2_cloud_init_script(),
                IamInstanceProfile = iam_instance_profile,
                BlockDeviceMappings=block_device_mappings,
                TagSpecifications=[
                    {
                        'ResourceType': 'instance',
                        'Tags': [{'Key': 'Name', 'Value': 'AWSEBSelfDestruct'}]
                    }
                ]
            )[0]
        except botocore.exceptions.ClientError as e:
            error_code = e.response['Error']['Code']
            if error_code == 'AccessDenied':
                print(f'Access denied! Please check your IAM permissions. \n   Error: {e}')
            else:
                print(f'Client Error: {e}')
            sys.exit(1)
        except Exception as e:
            print('Other Error: {e}')
            sys.exit(1)
    
        # Use a waiter to ensure the instance is running before trying to access its properties
        instance_id = instance.id    

        # tag the instance for cost explorer 
        tag = {
            'Key': 'INSTANCE_ID',
            'Value': instance_id
        }
        try:
            ec2.create_tags(Resources=[instance_id], Tags=[tag])
        except Exception as e:
            self.cfg.printdbg('Error creating Tags: {e}')
            
        print(f'Launching instance {instance_id} ... please wait ...')    
        
        max_wait_time = 300  # seconds
        delay_time = 10  # check every 10 seconds, adjust as needed
        max_attempts = max_wait_time // delay_time

        waiter = client.get_waiter('instance_running')
        progress = self._create_progress_bar(max_attempts)

        for attempt in range(max_attempts):
            try:
                waiter.wait(InstanceIds=[instance_id], WaiterConfig={'Delay': delay_time, 'MaxAttempts': 1})
                progress(attempt)
                break
            except botocore.exceptions.WaiterError:
                progress(attempt)
                continue
        print('')
        instance.reload()        

        grpid = self._ec2_create_and_attach_security_group(instance_id, profile)
        if grpid:
            print(f'Security Group "{grpid}" attached.') 
        else:
            print('No Security Group ID created.')
        instance.wait_until_running()
        print(f'Instance IP: {instance.public_ip_address}')

        self.cfg.write('cloud', 'ec2_last_instance', instance.public_ip_address)

        return instance_id, instance.public_ip_address

    def ec2_terminate_instance(self, ip, profile=None):
        # terminate instance  
        # with ephemeral (local) disk for a temporary restore 

        session = boto3.Session(profile_name=profile) if profile else boto3.Session()        
        ec2 = session.client('ec2')
        #ips = self.ec2_list_ips(self, 'Name', 'AWSEBSelfDestruct')    
        # Use describe_instances with a filter for the public IP address to find the instance ID
        filters = [{
            'Name': 'network-interface.addresses.association.public-ip',
            'Values': [ip]
        }]

        if not ip.startswith('i-'): # this an ip and not an instance ID
            try:
                response = ec2.describe_instances(Filters=filters)        
            except botocore.exceptions.ClientError as e: 
                print(f'Error: {e}')
                return False
            # Check if any instances match the criteria
            instances = [instance for reservation in response['Reservations'] for instance in reservation['Instances']]        
            if not instances:
                print(f"No EC2 instance found with public IP: {ip}")
                return 
            # Extract instance ID from the instance
            instance_id = instances[0]['InstanceId']
        else:    
            instance_id = ip
        # Terminate the instance
        ec2.terminate_instances(InstanceIds=[instance_id])
        
        print(f"EC2 Instance {instance_id} ({ip}) is being terminated !")

    def ec2_list_instances(self, tag_name, tag_value, profile=None):
        """
        List all IP addresses of running EC2 instances with a specific tag name and value.
        :param tag_name: The name of the tag
        :param tag_value: The value of the tag
        :return: List of IP addresses
        """
        session = boto3.Session(profile_name=profile) if profile else boto3.Session()
        ec2 = session.client('ec2')        
        
        # Define the filter
        filters = [
            {
                'Name': 'tag:' + tag_name,
                'Values': [tag_value]
            },
            {
                'Name': 'instance-state-name',
                'Values': ['running']
            }
        ]
        
        # Make the describe instances call
        try:
            response = ec2.describe_instances(Filters=filters)        
        except botocore.exceptions.ClientError as e:
            error_code = e.response['Error']['Code']
            if error_code == 'AccessDenied':
                self.cfg.printdbg(f'Access denied! Please check your IAM permissions. \n   Error: {e}')
            else:
                print(f'Client Error: {e}')
            return []            
        #An error occurred (AuthFailure) when calling the DescribeInstances operation: AWS was not able to validate the provided access credentials
        ilist = []    
        # Extract IP addresses
        for reservation in response['Reservations']:
            for instance in reservation['Instances']:
                row = [instance['PublicIpAddress'], 
                       instance['InstanceId'], 
                       instance['InstanceType']]
                ilist.append(row)                
        return ilist

    def ssh_execute(self, user, host, command=None):
        """Execute an SSH command on the remote server."""
        SSH_OPTIONS = "-o StrictHostKeyChecking=no"
        key_path = os.path.join(self.cfg.config_root,'cloud',f'{self.cfg.ssh_key_name}.pem')
        cmd = f"ssh {SSH_OPTIONS} -i '{key_path}' {user}@{host}"
        if command:
            cmd += f" '{command}'"
            try:
                result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
                return result
            except:
                print(f'Error executing "{cmd}."')
        else:
            subprocess.run(cmd, shell=True, capture_output=False, text=True)
        self.cfg.printdbg(f'ssh command line: {cmd}')
        return None
                
    def ssh_upload(self, user, host, local_path, remote_path, is_string=False):
        """Upload a file to the remote server using SCP."""
        SSH_OPTIONS = "-o StrictHostKeyChecking=no"
        key_path = os.path.join(self.cfg.config_root,'cloud',f'{self.cfg.ssh_key_name}.pem')
        if is_string:
            # the local_path is actually a string that needs to go into temp file 
            with tempfile.NamedTemporaryFile(mode='w+', delete=False) as temp:
                temp.write(local_path)
                local_path = temp.name
        cmd = f"scp {SSH_OPTIONS} -i '{key_path}' {local_path} {user}@{host}:{remote_path}"        
        try:
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            if is_string:
                os.remove(local_path)
            return result       
        except:
            print(f'Error executing "{cmd}."')
        return None

    def ssh_download(self, user, host, remote_path, local_path):
        """Upload a file to the remote server using SCP."""
        SSH_OPTIONS = "-o StrictHostKeyChecking=no"
        key_path = os.path.join(self.cfg.config_root,'cloud',f'{self.cfg.ssh_key_name}.pem')
        cmd = f"scp {SSH_OPTIONS} -i '{key_path}' {user}@{host}:{remote_path} {local_path}"        
        try:
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            return result
        except:
            print(f'Error executing "{cmd}."')
        return None            
    
    def send_email_ses(self, sender, to, subject, body, profile=None):
        # Using AWS ses service to send emails
        session = boto3.Session(profile_name=profile) if profile else boto3.Session()
        ses = session.client("ses")

        ses_verify_requests_sent = []
        if not sender:
            sender = self.cfg.read('general', 'email')
        if not to:
            to = self.cfg.read('general', 'email')
        if not to or not sender:
            print('from and to email addresses cannot be empty')
            return False
        ret = self.cfg.read('cloud', 'ses_verify_requests_sent')
        if isinstance(ret, list):
            ses_verify_requests_sent = ret
        else:
            ses_verify_requests_sent.append(ret)

        verified_email_addr = []
        try:
            response = ses.list_verified_email_addresses()
            verified_email_addr = response.get('VerifiedEmailAddresses', [])
        except botocore.exceptions.ClientError as e:
            error_code = e.response['Error']['Code']
            if error_code == 'AccessDenied':
                self.cfg.printdbg(f'Access denied to SES advanced features! Please check your IAM permissions. \nError: {e}')
            else:
                print(f'Client Error: {e}')
        except Exception as e:
            print(f'Other Error: {e}')
    
        checks = [sender, to]
        checks = list(set(checks)) # remove duplicates
        checked = []

        try:
            for check in checks:
                if check not in verified_email_addr and check not in ses_verify_requests_sent:
                    response = ses.verify_email_identity(EmailAddress=check)
                    checked.append(check)
                    print(f'{check} was used for the first time, verification email sent.')
                    print('Please have {check} check inbox and confirm email from AWS.\n')

        except botocore.exceptions.ClientError as e:
            error_code = e.response['Error']['Code']
            if error_code == 'AccessDenied':
                self.cfg.printdbg(f'Access denied to SES advanced features! Please check your IAM permissions. \nError: {e}')
            else:
                print(f'Client Error: {e}')
        except Exception as e:
            print(f'Other Error: {e}')
        
        self.cfg.write('cloud', 'ses_verify_requests_sent', checked)
        try:
            response = ses.send_email(
                Source=sender,
                Destination={
                    'ToAddresses': [to]
                },
                Message={
                    'Subject': {
                        'Data': subject
                    },
                    'Body': {
                        'Text': {
                            'Data': body
                        }
                    }
                }
            )
            print(f'Sent email "{subject}" to {to}!')
        except botocore.exceptions.ClientError as e:
            error_code = e.response['Error']['Code']
            if error_code == 'MessageRejected':
                print(f'Message was rejected, Error: {e}')
            elif error_code == 'AccessDenied':
                self.cfg.printdbg(f'Access denied to SES advanced features! Please check your IAM permissions. \nError: {e}')
                if not args.debug:
                    print (' Cannot use SES email features to send you status messages: AccessDenied')                
            else:
                print(f'Client Error: {e}')
            return False
        except Exception as e:
            print(f'Other Error: {e}')
            return False
        return True

        # The below AIM policy is needed if you do not want to confirm 
        # each and every email you want to send to. 

        # iam = boto3.client('iam')
        # policy_document = {
        #     "Version": "2012-10-17",
        #     "Statement": [
        #         {
        #             "Effect": "Allow",
        #             "Action": [
        #                 "ses:SendEmail",
        #                 "ses:SendRawEmail"
        #             ],
        #             "Resource": "*"
        #         }
        #     ]
        # }

        # policy_name = 'SES_SendEmail_Policy'

        # policy_arn = self._ec2_create_or_get_iam_policy(
        #     policy_name, policy_document, profile)

        # username = 'your_iam_username'  # Change this to the username you wish to attach the policy to

        # response = iam.attach_user_policy(
        #     UserName=username,
        #     PolicyArn=policy_arn
        # )

        # print(f"Policy {policy_arn} attached to user {username}")

    def send_ec2_costs(self, instance_id, profile=None):
        pass


    def _ec2_create_iam_costexplorer_ses(self, instance_id ,profile=None):
        session = boto3.Session(profile_name=profile) if profile else boto3.Session()
        iam = session.client('iam')
        ec2 = session.client('ec2')

        # Define the policy
        policy_document = {
            "Version": "2012-10-17",
            "Statement": [
                {
                    "Effect": "Allow",
                    "Action": "ses:SendEmail",
                    "Resource": "*"
                },                
                {
                    "Effect": "Allow",
                    "Action": [
                        "ce:*",              # Permissions for Cost Explorer
                        "ce:GetCostAndUsage", # all
                        "ec2:Describe*",     # Basic EC2 read permissions
                    ],
                    "Resource": "*"
                }
            ]
        }

        # Step 1: Create the policy in IAM
        policy_name = "CostExplorerSESPolicy"

        policy_arn = self._ec2_create_or_get_iam_policy(
            policy_name, policy_document, profile)
                
        # Step 2: Retrieve the IAM instance profile attached to the EC2 instance
        response = ec2.describe_instances(InstanceIds=[instance_id])
        instance_data = response['Reservations'][0]['Instances'][0]
        if 'IamInstanceProfile' not in instance_data:
            print(f"No IAM Instance Profile attached to the instance: {instance_id}")
            return False

        instance_profile_arn = response['Reservations'][0]['Instances'][0]['IamInstanceProfile']['Arn']

        # Extract the instance profile name from its ARN
        instance_profile_name = instance_profile_arn.split('/')[-1]

        # Step 3: Fetch the role name from the instance profile
        response = iam.get_instance_profile(InstanceProfileName=instance_profile_name)
        role_name = response['InstanceProfile']['Roles'][0]['RoleName']

        # Step 4: Attach the desired policy to the role
        try:
            iam.attach_role_policy(
                RoleName=role_name,
                PolicyArn=policy_arn
            )
            print(f"Policy {policy_arn} attached to role {role_name}")
        except iam.exceptions.NoSuchEntityException:
            print(f"Role {role_name} does not exist!")
        except iam.exceptions.InvalidInputException as e:
            print(f"Invalid input: {e}")
        except Exception as e:
            print(f"Other Error: {e}")


    def _ec2_create_iam_self_destruct_role(self, profile):
        session = boto3.Session(profile_name=profile) if profile else boto3.Session()
        iam = session.client('iam')

        # Step 1: Create IAM policy
        policy_document = {
            "Version": "2012-10-17",
            "Statement": [
                {
                    "Effect": "Allow",
                    "Action": "ec2:TerminateInstances",
                    "Resource": "*",
                    "Condition": {
                        "StringEquals": {
                            "ec2:ResourceTag/Name": "AWSEBSelfDestruct"
                        }
                    }
                }
            ]
        }
        policy_name = 'SelfDestructPolicy'
        
        policy_arn = self._ec2_create_or_get_iam_policy(
            policy_name, policy_document, profile)

        # Step 2: Create an IAM role and attach the policy
        trust_relationship = {
            "Version": "2012-10-17",
            "Statement": [
                {
                    "Effect": "Allow",
                    "Principal": {
                        "Service": "ec2.amazonaws.com"
                    },
                    "Action": "sts:AssumeRole"
                }
            ]
        }

        role_name = 'SelfDestructRole'
        try:
            iam.create_role(
                RoleName=role_name,
                AssumeRolePolicyDocument=json.dumps(trust_relationship),
                Description='Allows EC2 instances to call AWS services on your behalf.'
            )
        except iam.exceptions.EntityAlreadyExistsException:
            print ('IAM SelfDestructRole already exists.')            

        iam.attach_role_policy(
            RoleName=role_name,
            PolicyArn=policy_arn
        )

        return True
    

    def _get_ec2_metadata(self, metadata_entry):

        # request 'local-hostname', 'public-hostname', 'local-ipv4', 'public-ipv4'

        # Define the base URL for the EC2 metadata service
        base_url = "http://169.254.169.254/latest/meta-data/"

        # Request a token with a TTL of 60 seconds
        token_url = "http://169.254.169.254/latest/api/token"
        token_headers = {"X-aws-ec2-metadata-token-ttl-seconds": "60"}
        try:
            token_response = requests.put(token_url, headers=token_headers, timeout=2)
        except Exception as e:
            print(f'Other Error: {e}')
            return ""
        token = token_response.text

        # Use the token to retrieve the specified metadata entry
        headers = {"X-aws-ec2-metadata-token": token}
        try:
            response = requests.get(base_url + metadata_entry, headers=headers, timeout=2)
        except Exception as e:
            print(f'Other Error: {e}')
            return ""            

        if response.status_code != 200:
            print(f"Error: Failed to retrieve metadata for entry: {metadata_entry}. HTTP Status Code: {response.status_code}")
            return ""
    
        return response.text
    
    def monitor_ec2(self):

        # if system is idle self-destroy 

        instance_id = self._get_ec2_metadata('instance-id')
        public_ip = self._get_ec2_metadata('public-ipv4')
        instance_type = self._get_ec2_metadata('instance-type')
        ami_id = self._get_ec2_metadata('ami-id')
        reservation_id = self._get_ec2_metadata('reservation-id')
        
        nowstr = datetime.datetime.now().strftime('%H:%M:%S')
        print(f'aws-eb-monitor ({nowstr}): {public_ip} ({instance_id}, {instance_type}, {ami_id}, {reservation_id}) ... ')

        if self._monitor_is_idle():
            # This machine was idle for a long time, destroy it
            print(f'aws-eb-monitor ({nowstr}): Destroying current idling machine {public_ip} ({instance_id}) ...')
            if public_ip:
                body_text = "Instance was detected as idle and terminated"
                self.send_email_ses("", "", f'Terminating idle instance {public_ip} ({instance_id})', body_text)
                self.ec2_terminate_instance(public_ip)
                return True 
            else:
                print('Could not retrieve metadata (IP)')
                return False 
            
        current_time = datetime.datetime.now().time()
        start_time = datetime.datetime.strptime("23:00:00", "%H:%M:%S").time()
        end_time = datetime.datetime.strptime("23:59:59", "%H:%M:%S").time()    
        if start_time >= current_time or current_time > end_time:
            # only run cost emails once a day 
            return True 

        monthly_cost, monthly_unit, daily_costs_by_instance, user_monthly_cost, user_monthly_unit, \
            user_daily_cost, user_daily_unit, user_name = self._monitor_get_ec2_costs()
        
        body = []
        body.append(f"{monthly_cost:.2f} {monthly_unit} total account cost for the current month.")
        body.append(f"{user_monthly_cost:.2f} {user_monthly_unit} cost of user {user_name} for the current month.")
        body.append(f"{user_daily_cost:.2f} {user_daily_unit} cost of user {user_name} in the last 24 hours.")
        body.append("Cost for each EC2 instance type in the last 24 hours:")
        for instance_t, (cost, unit) in daily_costs_by_instance.items():
            if instance_t != 'NoInstanceType':
                body.append(f"  {instance_t:12}: ${cost:.2f} {unit}")
        body_text = "\n".join(body)
        self.send_email_ses("", "", f'AWS-EB AWS cost report ({instance_id})', body_text)

    def _monitor_users_logged_in(self):
        """Check if any users are logged in."""
        try:
            output = subprocess.check_output(['who']).decode('utf-8')
            if output:
                print('aws-eb-monitor: Not idle, logged in:', output)
                return True  # Users are logged in
            return False
        except Exception as e:
            print(f'Other Error: {e}')
            return True
        
    def _monitor_is_idle(self, interval=60, min_idle_cnt=72):

        # each run checks idle time for 60 seconds (interval)
        # if the system has been idle for 72 consecutive runs
        # the fucntion will return idle state after 3 days 
        # if the cron job is running hourly 

        # Constants
        CPU_THRESHOLD = 20  # percent
        NET_READ_THRESHOLD = 1000  # bytes per second
        NET_WRITE_THRESHOLD = 1000  # bytes per second
        DISK_WRITE_THRESHOLD = 100000  # bytes per second
        PROCESS_CPU_THRESHOLD = 10  # percent (for individual processes)
        PROCESS_MEM_THRESHOLD = 10  # percent (for individual processes)
        DISK_WRITE_EXCLUSIONS = ["systemd", "systemd-journald", \
                                "chronyd", "sshd", "auditd" , "agetty"]
    
        # Not idle if any users are logged in 
        if self._monitor_users_logged_in():
            print(f'aws-eb-monitor: Not idle: user(s) logged in')
            #return self._monitor_save_idle_state(False, min_idle_cnt)        
        
        # CPU, Time I/O and Network Activity 
        io_start = psutil.disk_io_counters()
        net_start = psutil.net_io_counters()
        cpu_percent = psutil.cpu_percent(interval=interval)
        io_end = psutil.disk_io_counters()
        net_end = psutil.net_io_counters()

        print(f'aws-eb-monitor: Current CPU% {cpu_percent}')

        # Check CPU Utilization
        if cpu_percent > CPU_THRESHOLD:
            print(f'aws-eb-monitor: Not idle: CPU% {cpu_percent}')
            #return self._monitor_save_idle_state(False, min_idle_cnt)

        # Check I/O Activity
        write_diff = io_end.write_bytes - io_start.write_bytes
        write_per_second = write_diff / interval

        if write_per_second > DISK_WRITE_THRESHOLD:
            for proc in psutil.process_iter(['name']):
                if proc.info['name'] in DISK_WRITE_EXCLUSIONS:
                    continue
                try:
                    if proc.io_counters().write_bytes > 0:
                        print(f'aws-eb-monitor:io bytes written: {proc.io_counters().write_bytes}')
                        return self._monitor_save_idle_state(False, min_idle_cnt)
                except (psutil.NoSuchProcess, psutil.AccessDenied):
                    pass

        # Check Network Activity
        bytes_sent_diff = net_end.bytes_sent - net_start.bytes_sent
        bytes_recv_diff = net_end.bytes_recv - net_start.bytes_recv

        bytes_sent_per_second = bytes_sent_diff / interval
        bytes_recv_per_second = bytes_recv_diff / interval

        if bytes_sent_per_second > NET_WRITE_THRESHOLD or \
                            bytes_recv_per_second > NET_READ_THRESHOLD:
            print(f'aws-eb-monitor:net bytes recv: {bytes_recv_per_second}')
            return self._monitor_save_idle_state(False, min_idle_cnt)
            
        # Examine Running Processes for CPU and Memory Usage
        for proc in psutil.process_iter(['name', 'cpu_percent', 'memory_percent']):
            if proc.info['name'] not in DISK_WRITE_EXCLUSIONS:
                if proc.info['cpu_percent'] > PROCESS_CPU_THRESHOLD:
                    print(f'aws-eb-monitor: Not idle: CPU% {proc.info["cpu_percent"]}')
                    #return False
                # disabled this idle checker 
                #if proc.info['memory_percent'] > PROCESS_MEM_THRESHOLD:
                #    print(f'aws-eb-monitor: Not idle: MEM% {proc.info["memory_percent"]}')
                #    return False

        # Write idle state and read consecutive idle hours
        print(f'aws-eb-monitor: Idle state detected')
        return self._monitor_save_idle_state(True, min_idle_cnt)

    def _monitor_save_idle_state(self, is_system_idle, min_idle_cnt):
        IDLE_STATE_FILE = os.path.join(os.getenv('TMPDIR', '/tmp'), 
                            'aws-eb_idle_state.txt')
        with open(IDLE_STATE_FILE, 'a') as file:
            file.write('1\n' if is_system_idle else '0\n')        
        with open(IDLE_STATE_FILE, 'r') as file:
            states = file.readlines()        
        count = 0
        for state in reversed(states):
            if state.strip() == '1':
                count += 1
            else:
                break
        return count >= min_idle_cnt        

    def _monitor_get_ec2_costs(self, profile=None):
        session = boto3.Session(profile_name=profile) if profile else boto3.Session()
        
        # Set up boto3 client for Cost Explorer
        ce = session.client('ce')
        sts = session.client('sts')

        # Identify current user/account
        identity = sts.get_caller_identity()
        user_arn = identity['Arn']
        # Check if it's the root user
        is_root = ":root" in user_arn

        # Dates for the current month and the last 24 hours
        today = datetime.datetime.today()
        first_day_of_month = datetime.datetime(today.year, today.month, 1).date()
        yesterday = (today - datetime.timedelta(days=1)).date()

        # Fetch EC2 cost of the current month
        monthly_response = ce.get_cost_and_usage(
            TimePeriod={
                'Start': str(first_day_of_month),
                'End': str(today.date())
            },
            Filter={
                'Dimensions': {
                    'Key': 'SERVICE',
                    'Values': ['Amazon Elastic Compute Cloud - Compute']
                }
            },
            Granularity='MONTHLY',
            Metrics=['UnblendedCost'],
        )
        monthly_cost = float(monthly_response['ResultsByTime'][0]['Total']['UnblendedCost']['Amount'])
        monthly_unit = monthly_response['ResultsByTime'][0]['Total']['UnblendedCost']['Unit']

        # If it's the root user, the whole account's costs are assumed to be caused by root.
        if is_root:
            user_name = 'root'
            user_monthly_cost = monthly_cost
            user_monthly_unit = monthly_unit
        else:
            # Assuming a tag `CreatedBy` (change as per your tagging system)
            user_name = user_arn.split('/')[-1]
            user_monthly_response = ce.get_cost_and_usage(
                TimePeriod={
                    'Start': str(first_day_of_month),
                    'End': str(today.date())
                },
                Filter={
                    "And": [
                        {
                            'Dimensions': {
                                'Key': 'SERVICE',
                                'Values': ['Amazon Elastic Compute Cloud - Compute']
                            }
                        },
                        {
                            'Tags': {
                                'Key': 'CreatedBy',
                                'Values': [user_name]
                            }
                        }
                    ]
                },
                Granularity='MONTHLY',
                Metrics=['UnblendedCost'],
            )
            user_monthly_cost = float(user_monthly_response['ResultsByTime'][0]['Total']['UnblendedCost']['Amount'])
            user_monthly_unit = user_monthly_response['ResultsByTime'][0]['Total']['UnblendedCost']['Unit']

        # Fetch cost of each EC2 instance type in the last 24 hours
        daily_response = ce.get_cost_and_usage(
            TimePeriod={
                'Start': str(yesterday),
                'End': str(today.date())
            },
            Filter={
                'Dimensions': {
                    'Key': 'SERVICE',
                    'Values': ['Amazon Elastic Compute Cloud - Compute']
                }
            },
            Granularity='DAILY',
            GroupBy=[{'Type': 'DIMENSION', 'Key': 'INSTANCE_TYPE'}],
            Metrics=['UnblendedCost'],
        )
        daily_costs_by_instance = {group['Keys'][0]: (float(group['Metrics']['UnblendedCost']['Amount']), group['Metrics']['UnblendedCost']['Unit']) for group in daily_response['ResultsByTime'][0]['Groups']}

        # Fetch cost caused by the current user in the last 24 hours
        if is_root:
            user_daily_cost = sum([cost[0] for cost in daily_costs_by_instance.values()])
            user_daily_unit = monthly_unit  # Using monthly unit since it should be the same for daily
        else:
            user_daily_response = ce.get_cost_and_usage(
                TimePeriod={
                    'Start': str(yesterday),
                    'End': str(today.date())
                },
                Filter={
                    "And": [
                        {
                            'Dimensions': {
                                'Key': 'SERVICE',
                                'Values': ['Amazon Elastic Compute Cloud - Compute']
                            }
                        },
                        {
                            'Tags': {
                                'Key': 'CreatedBy',
                                'Values': [user_name]
                            }
                        }
                    ]
                },
                Granularity='DAILY',
                Metrics=['UnblendedCost'],
            )
            user_daily_cost = float(user_daily_response['ResultsByTime'][0]['Total']['UnblendedCost']['Amount'])
            user_daily_unit = user_daily_response['ResultsByTime'][0]['Total']['UnblendedCost']['Unit']

        return monthly_cost, monthly_unit, daily_costs_by_instance, user_monthly_cost, \
               user_monthly_unit, user_daily_cost, user_daily_unit, user_name
    
class ConfigManager:
    # we write all config entries as files to '~/.config'
    # to make it easier for bash users to read entries 
    # with a simple var=$(cat ~/.config/aws-eb/section/entry)
    # entries can be strings, lists that are written as 
    # multi-line files and dictionaries which are written to json

    def __init__(self, args):
        self.args = args
        self.home_dir = os.path.expanduser('~')
        self.config_root_local = os.path.join(self.home_dir, '.config', 'aws-eb')
        self.config_root = self._get_config_root()
        self.binfolder = self.read('general', 'binfolder').replace(self.home_dir,'~')
        self.binfolderx = os.path.expanduser(self.binfolder)
        self.homepaths = self._get_home_paths()
        self.awscredsfile = os.path.join(self.home_dir, '.aws', 'credentials')
        self.awsconfigfile = os.path.join(self.home_dir, '.aws', 'config')
        self.awsconfigfileshr = os.path.join(self.config_root, 'aws_config')
        self.bucket = self.read('general','bucket')
        self.archiveroot = self.read('general','archiveroot')
        self.archivepath = f'{self.bucket}/{self.archiveroot}'
        self.awsprofile = os.getenv('AWS_PROFILE', 'default')
        profs = self.get_aws_profiles()
        if "aws" in profs:
            self.awsprofile = os.getenv('AWS_PROFILE', 'aws')
        elif "AWS" in profs:
            self.awsprofile = os.getenv('AWS_PROFILE', 'AWS')
        if hasattr(self.args, "awsprofile") and args.awsprofile:
            self.awsprofile = self.args.awsprofile
        self.aws_region = self.get_aws_region(self.awsprofile)
        self.envrn = os.environ.copy()
        if not self._set_env_vars(self.awsprofile):
            self.awsprofile = ''
        self.ssh_key_name = 'aws-eb-ec2'
        
    def _set_env_vars(self, profile):
        
        # Read the credentials file
        config = configparser.ConfigParser()
        config.read(self.awscredsfile)
        self.aws_region = self.get_aws_region(profile)

        if not config.has_section(profile):
            if self.args.debug:
                print (f'~/.aws/credentials has no section for profile {profile}')
            return False
        if not config.has_option(profile, 'aws_access_key_id'):
            if self.args.debug:
                print (f'~/.aws/credentials has no entry aws_access_key_id in section/profile {profile}')
            return False
        
        # Get the AWS access key and secret key from the specified profile
        aws_access_key_id = config.get(profile, 'aws_access_key_id')
        aws_secret_access_key = config.get(profile, 'aws_secret_access_key')

        # Set the environment variables for creds
        os.environ['AWS_ACCESS_KEY_ID'] = aws_access_key_id
        os.environ['AWS_SECRET_ACCESS_KEY'] = aws_secret_access_key
        os.environ['AWS_PROFILE'] = profile
        self.envrn['AWS_ACCESS_KEY_ID'] = aws_access_key_id
        self.envrn['AWS_SECRET_ACCESS_KEY'] = aws_secret_access_key
        self.envrn['AWS_PROFILE'] = profile
        self.envrn['RCLONE_S3_ACCESS_KEY_ID'] = aws_access_key_id
        self.envrn['RCLONE_S3_SECRET_ACCESS_KEY'] = aws_secret_access_key

        if profile in ['default', 'AWS', 'aws']:
            # Set the environment variables for AWS 
            self.envrn['RCLONE_S3_PROVIDER'] = 'AWS'
            self.envrn['RCLONE_S3_REGION'] = self.aws_region
            self.envrn['RCLONE_S3_LOCATION_CONSTRAINT'] = self.aws_region
            self.envrn['RCLONE_S3_STORAGE_CLASS'] = self.read('general','s3_storage_class')
            os.environ['RCLONE_S3_STORAGE_CLASS'] = self.read('general','s3_storage_class')
        else:
            prf=self.read('profiles',profile)
            self.envrn['RCLONE_S3_ENV_AUTH'] = 'true'
            self.envrn['RCLONE_S3_PROFILE'] = profile
            if isinstance(prf,dict):  # profile={'name': '', 'provider': '', 'storage_class': ''}
                self.envrn['RCLONE_S3_PROVIDER'] = prf['provider']
                self.envrn['RCLONE_S3_ENDPOINT'] = self._get_aws_s3_session_endpoint_url(profile)
                self.envrn['RCLONE_S3_REGION'] = self.aws_region
                self.envrn['RCLONE_S3_LOCATION_CONSTRAINT'] = self.aws_region
                self.envrn['RCLONE_S3_STORAGE_CLASS'] = prf['storage_class']
                os.environ['RCLONE_S3_STORAGE_CLASS'] = prf['storage_class']

        return True

    def sversion(self, version_str):
        """
        Parse a semantic versioning string into a tuple of integers.
        Args:
        version_str (str): A string representing the version, e.g., "8.2.45".
        Returns:
        tuple: A tuple of integers representing the major, minor, and patch versions.
        """
        parts = version_str.split('.')
        version = []
        for part in parts:
            if part.isdigit():
                version.append(int(part))  # Convert numeric strings to integers
            else:
                version.append(part)  # Keep non-numeric strings as is
        return tuple(version)
    
    def get_os_release_info(self):
        try:
            # Initialize the values
            os_id = ""
            version_id = ""
            # Open the file and read line by line
            with open('/etc/os-release', 'r') as f:
                for line in f:               
                    # Split the line into key and value
                    if line.startswith("ID="):                        
                        os_id = line.strip().split('=')[1].strip('"')
                    elif line.startswith("VERSION_ID="):                    
                        version_id = line.strip().split('=')[1].strip('"')
            return os_id, version_id
        except Exception as e:
            # Return two empty strings in case of an error
            return "", ""

    def _get_home_paths(self):
        path_dirs = os.environ['PATH'].split(os.pathsep)
        # Filter the directories in the PATH that are inside the home directory
        dirs_inside_home = {
            directory for directory in path_dirs
            if directory.startswith(self.home_dir) and os.path.isdir(directory)
        }
        return sorted(dirs_inside_home, key=len)  

    def _get_config_root(self):
        theroot=self.config_root_local
        rootfile = os.path.join(theroot, 'config_root')
        if os.path.exists(rootfile):
            with open(rootfile, 'r') as myfile:
                theroot = myfile.read().strip()
                if not os.path.isdir(theroot):
                    if not self.ask_yes_no(f'{rootfile} points to a shared config that does not exist. Do you want to configure {theroot} now?'):
                        print (f"Please remove file {rootfile} to continue with a single user config.")
                        sys.exit(1)
                        #raise FileNotFoundError(f'Config root folder "{theroot}" not found. Please remove {rootfile}')
        return theroot

    def _get_section_path(self, section):
        return os.path.join(self.config_root, section)

    def _get_entry_path(self, section, entry):
        if section:
            section_path = self._get_section_path(section)
            return os.path.join(section_path, entry)
        else:
            return os.path.join(self.config_root, entry)

    def replace_symlinks_with_realpaths(self, folders):
        cleaned_folders = []
        for folder in folders:
            try:
                # Split the path into its components
                folder = os.path.expanduser(folder)
                #print('expanduser folder:', folder)
                #print('real path:', os.path.realpath(folder))
                cleaned_folders.append(os.path.realpath(folder))
            except Exception as e:
                print(f"Error processing '{folder}': {e}")           
        self.printdbg('cleaned_folders:', cleaned_folders)
        return cleaned_folders            
        
    def printdbg(self, *args, **kwargs):
        # use inspect to get the name of the calling function
        if self.args.debug:
            current_frame = inspect.currentframe()
            calling_function = current_frame.f_back.f_code.co_name 
            print(f' DBG {calling_function}():', args, kwargs)

    def prompt(self, question, defaults=None, type_check=None):
        # Prompts for user input and writes it to config. 
        # defaults are up to 3 pipe separated strings: 
        # if there is only one string this is the default 
        #
        # if there are 2 strings they represent section 
        # and key name of the config entry and if there are 
        # 3 strings the last 2 represent section 
        # and key name of the config file and the first is
        # the default if section and key name are empty
        #
        # if defaults is a python list it will assign a number
        # to each list element and prompt the user for one
        # of the options
        default=''
        section=''
        key=''
        if not question.endswith(':'):
            question += ':'
        question = f"*** {question} ***"
        if isinstance(defaults, list):
            print(question)
            for i, option in enumerate(defaults, 1):
                print(f'  ({i}) {option}')           
            while True:
                selected = input("  Enter the number of your selection: ")
                if selected.isdigit() and 1 <= int(selected) <= len(defaults):
                    return defaults[int(selected) - 1]
                else:
                    print("  Invalid selection. Please enter a number from the list.")
        elif defaults is not None:
            deflist=defaults.split('|')
            if len(deflist) == 3:
                section=deflist[1]
                key=deflist[2]
                default = self.read(section, key)
                if not default:
                    default = deflist[0]
            elif len(deflist) == 2:
                section=deflist[0]
                key=deflist[1]
                default = self.read(section, key)                                
            elif len(deflist) == 1:
                default = deflist[0]
            #if default:
            question += f"\n  [Default: {default}]"
        else:
            question += f"\n  [Default: '']"
        while True:
            #user_input = input(f"\033[93m{question}\033[0m ")
            user_input = input(f"{question} ")
            if not user_input:
                if default is not None:
                    if section:
                        self.write(section,key,default)                    
                    return default
                else:
                    print("Please enter a value.")
            else:
                if type_check == 'number':
                    try:
                        if '.' in user_input:
                            value = float(user_input)
                        else:
                            value = int(user_input)
                        if section:
                            self.write(section,key,value)
                        return value
                    except ValueError:
                        print("Invalid input. Please enter a number.")
                elif type_check == 'string':
                    if not user_input.isnumeric():
                        if section:
                            self.write(section,key,user_input)
                        return user_input
                    else:
                        print("Invalid input. Please enter a string not a number")
                else:
                    if section:
                        self.write(section,key,user_input)
                    return user_input

    def ask_yes_no(self, question, default="yes"):
        valid = {"yes": True, "y": True, "no": False, "n": False}

        if default is None:
            prompt = " [y/n] "
        elif default == "yes":
            prompt = " [Y/n] "
        elif default == "no":
            prompt = " [y/N] "
        else:
            raise ValueError("invalid default answer: '%s'" % default)

        while True:
            print(question + prompt, end="")
            choice = input().lower()
            if default and not choice:
                return valid[default]
            elif choice in valid:
                return valid[choice]
            else:
                print("Please respond with 'yes' or 'no' (or 'y' or 'n').")


    def add_cron_job(self, cmd, minute, hour='*', day_of_month='*', month='*', day_of_week='*'):
        # CURRENTLY INACTIVE
        if not minute:
            print('You must set the minute (1-60) explicily')
            return False 
        with tempfile.NamedTemporaryFile(delete=False) as temp:
            # Dump the current crontab to the temporary file
            try:
                os.system('crontab -l > {}'.format(temp.name))
            except Exception as e:
                print(f"Error: {e}")                

            # Add the new cron job to the temporary file
            cron_time = "{} {} {} {} {}".format(str(minute), hour, day_of_month, month, day_of_week)
            with open(temp.name, 'a') as file:
                file.write('{} {}\n'.format(cron_time, cmd))
            
            # Install the new crontab
            try:            
                os.system('crontab {}'.format(temp.name))
            except Exception as e:
                print(f"Error: {e}")                

            # Clean up by removing the temporary file
            os.unlink(temp.name)

        print("Cron job added!")

    def add_systemd_cron_job(self, cmd, minute, hour='*'):

        # Troubleshoot with: 
        #
        # journalctl -f --user-unit aws-eb-monitor.service
        # journalctl -f --user-unit aws-eb-monitor.timer
        # journalctl --since "5 minutes ago" | grep aws-eb-monitor

        SERVICE_CONTENT = textwrap.dedent(f"""
        [Unit]
        Description=Run AWS-EB-Monitor Cron Job

        [Service]
        Type=simple
        ExecStart={cmd}

        [Install]
        WantedBy=default.target
        """)

        TIMER_CONTENT = textwrap.dedent(f"""
        [Unit]
        Description=Run AWS-EB-Monitor Cron Job hourly

        [Timer]
        Persistent=true
        OnCalendar=*-*-* {hour}:{minute}:00
        #RandomizedDelaySec=300
        #FixedRandomDelay=true
        #OnBootSec=180
        #OnUnitActiveSec=3600
        Unit=aws-eb-monitor.service

        [Install]
        WantedBy=timers.target
        """)

        # Ensure the directory exists
        user_systemd_dir = os.path.expanduser("~/.config/systemd/user/")
        os.makedirs(user_systemd_dir, exist_ok=True)

        SERVICE_PATH = os.path.join(user_systemd_dir, "aws-eb-monitor.service")
        TIMER_PATH = os.path.join(user_systemd_dir, "aws-eb-monitor.timer")

        # Create service and timer files
        with open(SERVICE_PATH, "w") as service_file:
            service_file.write(SERVICE_CONTENT)

        with open(TIMER_PATH, "w") as timer_file:
            timer_file.write(TIMER_CONTENT)

        # Reload systemd and enable/start timer
        try:
            os.chdir(user_systemd_dir)
            os.system("systemctl --user daemon-reload")            
            os.system("systemctl --user enable aws-eb-monitor.service")
            os.system("systemctl --user enable aws-eb-monitor.timer")            
            os.system("systemctl --user start aws-eb-monitor.timer")            
            print("Systemd aws-eb-monitor.timer cron job started!")
        except Exception as e:
            print(f'Could not add systemd scheduler job, Error: {e}')

    def replicate_ini(self, section, src_file, dest_file):

        # copy an ini section from source to destination
        # sync values in dest that do not exist in src back to src
        # best used for sync of AWS profiles.
        # if section==ALL copy all but section called default

        if not os.path.exists(src_file):
            return

        # Create configparser objects
        src_parser = configparser.ConfigParser()
        dest_parser = configparser.ConfigParser()

        # Read source and destination files
        src_parser.read(src_file)
        dest_parser.read(dest_file)

        if section == 'ALL':
            sections = src_parser.sections()
            sections.remove('default') if 'default' in sections else None
        else:
            sections = [section]

        for section in sections:
            # Get the section from source and destination files
            src_section_data = dict(src_parser.items(section))
            dest_section_data = dict(dest_parser.items(section)) if dest_parser.has_section(section) else {}

            # If section does not exist in source or destination file, add it
            if not src_parser.has_section(section):
                src_parser.add_section(section)

            if not dest_parser.has_section(section):
                dest_parser.add_section(section)

            # Write the data into destination file
            for key, val in src_section_data.items():
                dest_parser.set(section, key, val)

            # Write the data into source file
            for key, val in dest_section_data.items():
                if key not in src_section_data:
                    src_parser.set(section, key, val)

        # Save the changes in the destination and source files
        with open(dest_file, 'w') as dest_configfile:
            dest_parser.write(dest_configfile)

        with open(src_file, 'w') as src_configfile:
            src_parser.write(src_configfile)

        if self.args.debug:
            print(f"Ini-section copied from {src_file} to {dest_file}")
            print(f"Missing entries in source from destination copied back to {src_file}")

    def get_aws_profiles(self):
        # get the full list of profiles from ~/.aws/ profile folder
        config = configparser.ConfigParser()        
        # Read the AWS config file ---- optional, we only require a creds file
        if os.path.exists(self.awsconfigfile):
            config.read(self.awsconfigfile)        
        # Read the AWS credentials file
        if os.path.exists(self.awscredsfile):
            config.read(self.awscredsfile)        
        # Get the list of profiles
        profiles = []
        for section in config.sections():
            profile_name = section.replace("profile ", "") #.replace("default", "default")
            profiles.append(profile_name)
        # convert list to set and back to list to remove dups
        return list(set(profiles))

    def create_aws_configs(self,access_key=None, secret_key=None, region=None):

        aws_dir = os.path.join(self.home_dir, ".aws")

        if not os.path.exists(aws_dir):
            os.makedirs(aws_dir)

        if not os.path.isfile(self.awsconfigfile):
            if region:
                print(f'\nAWS config file {self.awsconfigfile} does not exist, creating ...')            
                with open(self.awsconfigfile, "w") as config_file:
                    config_file.write("[default]\n")
                    config_file.write(f"region = {region}\n")
                    config_file.write("\n")
                    config_file.write("[profile aws]\n")
                    config_file.write(f"region = {region}\n")

        if not os.path.isfile(self.awscredsfile):
            print(f'\nAWS credentials file {self.awscredsfile} does not exist, creating ...')
            if not access_key: access_key = input("Enter your AWS access key ID: ")
            if not secret_key: secret_key = input("Enter your AWS secret access key: ")            
            with open(self.awscredsfile, "w") as credentials_file:
                credentials_file.write("[default]\n")
                credentials_file.write(f"aws_access_key_id = {access_key}\n")
                credentials_file.write(f"aws_secret_access_key = {secret_key}\n")
                credentials_file.write("\n")
                credentials_file.write("[aws]\n")
                credentials_file.write(f"aws_access_key_id = {access_key}\n")
                credentials_file.write(f"aws_secret_access_key = {secret_key}\n")
            os.chmod(self.awscredsfile, 0o600)

    def set_aws_config(self, profile, key, value, service=''):
        if key == 'endpoint_url': 
            if value.endswith('.amazonaws.com'):
                return False
            else:
                value = f'{value}\nsignature_version = s3v4'
        config = configparser.ConfigParser()
        config.read(os.path.expanduser("~/.aws/config"))
        section=profile
        if profile != 'default':
            section = f'profile {profile}'
        if not config.has_section(section):
            config.add_section(section)
        if service: 
            config.set(section, service, f"\n{key} = {value}\n")
        else:
            config.set(section, key, value)
        with open(os.path.expanduser("~/.aws/config"), 'w') as configfile:
            config.write(configfile)
        return True
    
    def get_aws_s3_endpoint_url(self, profile=None):
        # non boto3 method, use _get_aws_s3_session_endpoint_url instead
        if not profile:
            profile=self.awsprofile
        config = configparser.ConfigParser()
        config.read(os.path.expanduser('~/.aws/config'))
        prof = 'profile ' + profile
        if profile == 'default':
            prof = profile
        try:
            # We use the configparser's interpolation feature here to 
            # flatten the 's3' subsection into the 'profile test' section.
            s3_config_string = config.get(prof, 's3')
            s3_config = configparser.ConfigParser()
            s3_config.read_string("[s3_section]\n" + s3_config_string)
            endpoint_url = s3_config.get('s3_section', 'endpoint_url')
            return endpoint_url
        except (configparser.NoSectionError, configparser.NoOptionError):
            if self.args.debug:
                print("  No endpoint_url found in aws profile:", profile)
            return None
        
    def _get_aws_s3_session_endpoint_url(self, profile=None):
        # retrieve endpoint url through boto API, not configparser
        import botocore.session  # only botocore Session object has attribute 'full_config'
        if not profile:
            profile = self.awsprofile        
        session = botocore.session.Session(profile=profile) if profile else botocore.session.Session()        
        config = session.full_config
        s3_config = config["profiles"][profile].get("s3", {})
        endpoint_url = s3_config.get("endpoint_url", None)
        if self.args.debug:
            print('*** endpoint url ***:', endpoint_url)
        return endpoint_url

    def get_aws_region(self, profile=None):
        try:            
            session = boto3.Session(profile_name=profile) if profile else boto3.Session()
            if self.args.debug:
                print(f'* get_aws_region for profile {profile}:', session.region_name)
            return session.region_name
        except:
            if self.args.debug:
                print(f'  cannot retrieve AWS region for profile {profile}, no valid profile or credentials')
            return ""
            
    def get_domain_name(self):
        try:
            with open('/etc/resolv.conf', 'r') as file:
                content = file.readlines()
        except FileNotFoundError:
            return "mydomain.edu"
        tld = None
        for line in content:
            if line.startswith('search') or line.startswith('domain'):
                tokens = line.split()
                if len(tokens) > 1:
                    tld = tokens.pop()
                    break
        return tld if tld else "mydomain.edu"
    
    def get_time_zone(self):
        current_tz_str = 'America/Los_Angeles'
        try:        
            # Resolve the /etc/localtime symlink
            timezone_path = os.path.realpath("/etc/localtime")
            # Extract the time zone string by stripping off the prefix of the zoneinfo path
            current_tz_str = timezone_path.split("zoneinfo/")[-1]
            #import zoneinfo 
            #current_tz = zoneinfo.ZoneInfo(current_tz_str)\
            #current_time = datetime.datetime.now(current_tz)
        except Exception as e:
            print(f'Error: {e}')
            current_tz_str = 'America/Los_Angeles'
        return current_tz_str
        
    def write(self, section, entry, value):
        entry_path = self._get_entry_path(section, entry)
        os.makedirs(os.path.dirname(entry_path), exist_ok=True)
        if value == '""':
            os.remove(entry_path)
            return
        with open(entry_path, 'w') as entry_file:
            if isinstance(value, list):
                for item in value:
                    entry_file.write(f"{item}\n")
            elif isinstance(value, dict):
                json.dump(value, entry_file)
            else:
                entry_file.write(value)

    def read(self, section, entry):
        entry_path = self._get_entry_path(section, entry)
        if not os.path.exists(entry_path):
            return ""
            #raise FileNotFoundError(f'Config entry "{entry}" in section "{section}" not found.')
        with open(entry_path, 'r') as entry_file:
            try:
                return json.load(entry_file)                
            except json.JSONDecodeError:
                pass
            except:
                print('Error in ConfigManager.read()')
        with open(entry_path, 'r') as entry_file:
                content = entry_file.read().splitlines()
                if len(content) == 1:
                    return content[0].strip()
                else:
                    return content

    def delete(self, section, entry):
        entry_path = self._get_entry_path(section, entry)
        if not os.path.exists(entry_path):
            raise FileNotFoundError(f'Config entry "{entry}" in section "{section}" not found.')
        os.remove(entry_path)

    def delete_section(self, section):
        section_path = self._get_section_path(section)
        if not os.path.exists(section_path):
            raise FileNotFoundError(f'Config section "{section}" not found.')
        for entry in os.listdir(section_path):
            os.remove(os.path.join(section_path, entry))
        os.rmdir(section_path)

    def move_config(self,cfgfolder):
        if not cfgfolder and self.config_root == self.config_root_local:
                cfgfolder = self.prompt("Please enter the root where folder .config/aws-eb will be created.", 
                                    os.path.expanduser('~'))
        if cfgfolder:
            new_config_root = os.path.join(os.path.expanduser(cfgfolder),'.config','aws-eb')
        else:
            new_config_root = self.config_root
        old_config_root = self.config_root_local
        config_root_file = os.path.join(self.config_root_local,'config_root')
        
        if os.path.exists(config_root_file):
            with open(config_root_file, 'r') as f:
                old_config_root = f.read().strip()
        
        #print(old_config_root,new_config_root)
        if old_config_root == new_config_root:
            return True

        if not os.path.isdir(new_config_root):
            if os.path.isdir(old_config_root):
                shutil.move(old_config_root,new_config_root) 
                if os.path.isdir(old_config_root):
                    try:
                        os.rmdir(old_config_root)
                    except:
                        pass
                print(f'  AWS-EB config moved to "{new_config_root}"\n')
            os.makedirs(new_config_root,exist_ok=True)
            if os.path.exists(self.awsconfigfile):
                self.replicate_ini('ALL',self.awsconfigfile,os.path.join(new_config_root,'aws_config'))
                print(f'  ~/.aws/config replicated to "{new_config_root}/aws_config"\n')  

        self.config_root = new_config_root

        os.makedirs(old_config_root,exist_ok=True)
        with open(config_root_file, 'w') as f:
            f.write(self.config_root)
            print(f'  Switched configuration path to "{self.config_root}"')
        return True
    
    def wait_for_ssh_ready(self, hostname, port=22, timeout=60):
        start_time = time.time()
        while time.time() - start_time < timeout:
            s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            s.settimeout(3)  # Set a timeout on the socket operations
            result = s.connect_ex((hostname, port))
            if result == 0:
                s.close()
                return True
            else:
                time.sleep(5)  # Wait for 5 seconds before retrying
                s.close()
        print("Timeout reached without SSH server being ready.")
        return False

    def copy_compiled_binary_from_github(self,user,repo,compilecmd,binary,targetfolder):
        tarball_url = f"https://github.com/{user}/{repo}/archive/refs/heads/main.tar.gz"
        response = requests.get(tarball_url, stream=True, allow_redirects=True)
        response.raise_for_status()
        with tempfile.TemporaryDirectory() as tmpdirname:
            reposfolder=os.path.join(tmpdirname,  f"{repo}-main")
            with tarfile.open(fileobj=response.raw, mode="r|gz") as tar:
                tar.extractall(path=tmpdirname)
                reposfolder=os.path.join(tmpdirname,  f"{repo}-main")
                os.chdir(reposfolder)
                result = subprocess.run(compilecmd, shell=True)
                if result.returncode == 0:
                    print(f"Compilation successful: {compilecmd}")
                    shutil.copy2(binary, targetfolder, follow_symlinks=True)
                    if not os.path.exists(os.path.join(targetfolder, binary)):
                        print(f'Failed copying {binary} to {targetfolder}')                
                else:
                    print(f"Compilation failed: {compilecmd}")

    def copy_binary_from_zip_url(self,zipurl,binary,subwildcard,targetfolder):
        with tempfile.TemporaryDirectory() as tmpdirname:
            zip_file = os.path.join(tmpdirname,  "download.zip")
            response = requests.get(zipurl, verify=False, allow_redirects=True)
            with open(zip_file, 'wb') as f:
                f.write(response.content)
            with zipfile.ZipFile(zip_file, 'r') as zip_ref:
                zip_ref.extractall(tmpdirname)
            binpath = glob.glob(f'{tmpdirname}{subwildcard}{binary}')[0]
            shutil.copy2(binpath, targetfolder, follow_symlinks=True)
            if os.path.exists(os.path.join(targetfolder, binary)):
                os.chmod(os.path.join(targetfolder, binary), 0o775)
            else:    
                print(f'Failed copying {binary} to {targetfolder}')


# class SubcommandHelpFormatter(argparse.RawDescriptionHelpFormatter):
#     def _format_action(self, action):
#         parts = super(argparse.RawDescriptionHelpFormatter, self)._format_action(action)
#         if action.nargs == argparse.PARSER:
#             parts = "\n".join(parts.split("\n")[1:])
#         return parts


def parse_arguments():
    """
    Gather command-line arguments.
    """       
    parser = argparse.ArgumentParser(prog='aws-eb ',
        description='A (mostly) for building stuff on AWS  ' + \
                    'after finding folders in the file system that are worth archiving.')
    parser.add_argument( '--debug', '-d', dest='debug', action='store_true', default=False,
        help="verbose output for all commands")
    parser.add_argument('--profile', '-p', dest='awsprofile', action='store', default='', 
        help='which AWS profile in ~/.aws/ should be used. default="aws"')
    parser.add_argument('--version', '-v', dest='version', action='store_true', default=False, 
        help='print AWS-EB and Python version info')
    
    subparsers = parser.add_subparsers(dest="subcmd", help='sub-command help')

    # ***
    parser_config = subparsers.add_parser('config', aliases=['cnf'], 
        help=textwrap.dedent(f'''
            Bootstrap the configurtion, install dependencies and setup your environment.
            You will need to answer a few questions about your cloud and hpc setup.
        '''), formatter_class=argparse.RawTextHelpFormatter)
    parser_config.add_argument( '--monitor', '-m', dest='monitor', action='store', default='',
        metavar='<email@address.org>', help='setup aws-eb as a monitoring cronjob ' +
        'on an ec2 instance and notify an email address')

    # ***
    parser_launch = subparsers.add_parser('launch', aliases=['lau'],
        help=textwrap.dedent(f'''
            Launch EC2 instance, build new Easybuild packages and upload them to S3
        '''), formatter_class=argparse.RawTextHelpFormatter) 
    parser_launch.add_argument('--instance-type', '-i', dest='instancetype', action='store', default="",
        help='The EC2 instance type is auto-selected, but you can pick any other type here')    
    parser_launch.add_argument('--gpu-type', '-g', dest='gputype', action='store', default="",
        help='run --list to see available GPU types')       
    parser_launch.add_argument('--cpu-type', '-c', dest='cputype', action='store', default="",
        help='run --list to see available CPU types')        
    parser_launch.add_argument( '--list', '-l', dest='list', action='store_true', default=False,
        help="List CPU and GPU types")
    parser_launch.add_argument('--vcpus', '-v', dest='vcpus', type=int, action='store', default=4, 
        help='Number of cores to be allocated for the machine. (default=4)')    
    parser_launch.add_argument('--mem', '-m', dest='mem', type=int, action='store', default=8, 
        help='GB Memory allocated to instance  (default=8)')    
    parser_launch.add_argument( '--monitor', '-n', dest='monitor', action='store_true', default=False,
        help="Monitor EC2 server for cost and idle time.")
    parser_launch.add_argument( '--build', '-b', dest='build', action='store_true', default=False,
        help="Build the Easybuild packages on current system.")
    parser_launch.add_argument( '--bio-only', '-o', dest='bioonly', action='store_true', default=False,
        help="Build only life sciences (bio) packages and dependencies.")
    
    # ***
    parser_download = subparsers.add_parser('download', aliases=['dld'],
        help=textwrap.dedent(f'''
            Download built eb packages to /opt/eb
        '''), formatter_class=argparse.RawTextHelpFormatter) 
    parser_download.add_argument('--target', '-t', dest='target', action='store', default='/opt/eb', 
        metavar='<target_folder>', help='Download to other folder than default')    
    
    # ***
    parser_ssh = subparsers.add_parser('ssh', aliases=['scp'],
        help=textwrap.dedent(f'''
            Login to an AWS EC2 instance to which data was restored with the --ec2 option
        '''), formatter_class=argparse.RawTextHelpFormatter)
    parser_ssh.add_argument( '--list', '-l', dest='list', action='store_true', default=False,
        help="List running AWS-EB EC2 instances")        
    parser_ssh.add_argument('--terminate', '-t', dest='terminate', action='store', default='', 
        metavar='<hostname>', help='Terminate EC2 instance with this public IP Address.')    
    # parser_ssh.add_argument('--key', '-k', dest='key', action='store', default='', 
    #     help='pick a custom ssh key for your connection')
    parser_ssh.add_argument('sshargs', action='store', default=[], nargs='*',
        help='multiple arguments to ssh/scp such as hostname or user@hostname oder folder' +
               '')

    if len(sys.argv) == 1:
        parser.print_help(sys.stdout)               

    return parser.parse_args()

if __name__ == "__main__":
    if not sys.platform.startswith('linux'):
        print('This software currently only runs on Linux x64')
        sys.exit(1)
    try:
        args = parse_arguments()
        if main():
            sys.exit(0)
        else:
            sys.exit(1)
    except KeyboardInterrupt:
        print('\nExit !')
        try:
            sys.exit(0)
        except SystemExit:
            os._exit(0)
