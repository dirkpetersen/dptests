# AWS-HPC

AWS-HPC is a tool similar to Froster and SciBob. A command line tool that launches EC2 instances in the AWS spot market, installs scientific software, integrates S3 buckets as a FileSystem using JuiceFS. 
It also implements a subset of the sbatch/srun/salloc CLI from Slurm https://slurm.schedmd.com/. Unlike Slurm, it does not aim to spread jobs across multiple compute nodes as AWS has instances up to 224 cores (448 vcps) and very few users ever submit jobs that use more then 100 cores in a single job. Each submitted job will launch a new instance and work with a shared file system. 
It is implemented as a collaborative HPC system where userids and uid/gid are taken from Github to avoid naming conflicts. 





