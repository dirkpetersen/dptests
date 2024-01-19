
# Notes on building packages with EasyBuild and Slurm

Note: you may have to set Slurm input environment vars:
      https://slurm.schedmd.com/sbatch.html#lbAJ

* use the --job option to submit builds to the cluster 
eb myeasyconfig.eb --robot --job

* to build packages using the A30 GPU set environment
export SBATCH_GRES=gpu:a30:1

* some packages need explicit setting of OpenSSL
export OPENSSL_DIR=/app/eb/software/OpenSSL/1.1
 


