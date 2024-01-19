# .bash_profile

# Get the aliases and functions
if [ -f ~/.bashrc ]; then
	. ~/.bashrc
fi

# User specific environment and startup programs

# Slurm settings for software builds
export SBATCH_QOS=large
export SBATCH_MEM_PER_CPU=16G
export SBATCH_GRES=gpu:a30:1

# easybuild envionment for app builder, use 'eb --job' to submit build jobs to slurm  
export EASYBUILD_PREFIX=/app/eb
export EASYBUILD_ROBOT_PATHS=/app/eb/fh/fh_easyconfigs/:/app/eb/mcc/mcc_easyconfigs/ # Custom configs 
export EASYBUILD_ROBOT_PATHS=${EASYBUILD_ROBOT_PATHS}:/app/eb/devel/easybuild-easyconfigs/easybuild/easyconfigs # EB Develop branch
export EASYBUILD_ROBOT_PATHS=${EASYBUILD_ROBOT_PATHS}:~/.local/easybuild/easyconfigs # EB Production branch
export EASYBUILD_GITHUB_USER=scicomp-moffitt
export EASYBUILD_UPDATE_MODULES_TOOL_CACHE=True
export EASYBUILD_JOB_BACKEND=Slurm
export EASYBUILD_JOB_CORES=4
export EASYBUILD_JOB_MAX_WALLTIME=24
export EASYBUILD_PARALLEL=16
export EASYBUILD_DOWNLOAD_TIMEOUT=600
export EASYBUILD_JOB_OUTPUT_DIR=${EASYBUILD_PREFIX}/slurm-output
export EASYBUILD_BUILDPATH=/dev/shm/$USER
export EASYBUILD_CUDA_COMPUTE_CAPABILITIES="7.5,8.0" #A30 is 8.0, A40 is 8.6, H100 is 9.0
## Output is huge
#export EASYBUILD_LOGTOSTDOUT=True

# SSL hack needed for some EasyBuild Apps (e.g. for Rust)
export OPENSSL_DIR=/app/eb/software/OpenSSL/1.1

# other settings to manage /app
#export PATH=/app/lib/scadmin/bin:${PATH}

ml purge

#[ -z $SLURM_PTY_PORT ] && eval $(~/.local/bin/keychain --quiet --eval id_ed25519)
