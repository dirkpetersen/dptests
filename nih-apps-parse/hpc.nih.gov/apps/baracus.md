

document.querySelector('title').textContent = 'Baracus on Biowulf';
Baracus on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
[Swarm of jobs](#swarm) 
 |



BARACUS is a BIDS-compliant application that predicts brain age, based on combining data from cortical thickness, cortical surface area, and subcortical information



### Web site


* [Home page](https://github.com/BIDS-Apps/baracus)


Important Notes
* baracus on Biowulf is a singularity container built directly from docker://bids/baracus. However, users do not need to execute singularity directly or bind directories manually because we have provided a wrapper shell script that takes care of that. 
* Module Name: baracus (see [the modules page](/apps/modules.html) for more information)
* Executable Name: baracus* When using baracus we recommend the following:
	1. Limit I/O with central Biowulf directories (e.g., /data, /home) by making use of local disk (/lscratch). You can make use of lscratch by (a) copying the input dataset to /lscratch/$SLURM\_JOB\_ID and (b) by using the baracus flag out\_dir to write the ouput data to /lscratch/$SLURM\_JOB\_ID (remember to allocate enough space in /lscratch).
	* Limit the number of threads that baracus is allowed to use by including the --n\_cpus flag in the baracus command.
	* Profile/benchmark baracus jobs: We recommend making sure a given baracus job can scale before launching a large number of jobs. You can do this by profiling/benchmarking small jobs (e.g., a swarm of 3 baracus commands), then monitoring the jobs by using either the [user dashboard](https://hpc.nih.gov/dashboard/) or the commands jobhist, sjobs, squeue, and jobload (see [biowulf utilities](https://hpc.nih.gov/docs/biowulf_tools.html)). There are resources prepared by the HPC staff that go over how to monitor your jobs for benchmarking and profiling ([video](https://youtu.be/fLMJ8-t5bm4), and [slides](https://hpc.nih.gov/training/handouts/Effective_batch_system.pdf)). Once you have profiled a swarm with a few jobs and you have refined the memory, CPU (and walltime) requirements, it is probably safe to (gradually) increase the number of baracus jobs. For many scientific pipelines one has no way of knowing in advance how much memory or CPUs will be actually required in an HPC environment. This is why it is very important to profile/benchmark.



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load baracus**

[user@cn3144 ~]$ **baracus -h**
usage: run_brain_age_bids.py [-h]
                             [--participant_label PARTICIPANT_LABEL [PARTICIPANT_LABEL ...]]
                             [--freesurfer_dir FREESURFER_DIR]
                             [--models {Liem2016__OCI_norm,Liem2016__full_2samp_training} [{Liem2016__OCI_norm,Liem2016__full_2samp_training} ...]]
                             [--skip_missing] --license_key LICENSE_KEY
                             [--n_cpus N_CPUS] [-v]
                             bids_dir out_dir {participant,group}

BARACUS: Brain-Age Regression Analysis and Computation Utility Software. BIDS
mode. You specify a BIDS-formatted freesurfer folder as input. All data is
extracted automatiacally from that folder.

positional arguments:
  bids_dir              The directory with the input dataset formatted
                        according to the BIDS standard.
  out_dir               Results are put into {out_dir}/baracus.
  {participant,group}   Level of the analysis that will be performed.
                        "participant": predicts single subject brain age,
                        "group": collects single subject predictions.

optional arguments:
  -h, --help            show this help message and exit
  --participant_label PARTICIPANT_LABEL [PARTICIPANT_LABEL ...]
                        The label of the participant that should be analyzed.
                        The label corresponds to sub- from
 the BIDS spec (so it does not include "sub-"). If this
 parameter is not provided all subjects should be
 analyzed. Multiple participants can be specified with
 a space separated list.
 --freesurfer\_dir FREESURFER\_DIR
 Folder with FreeSurfer subjects formatted according to
 BIDS standard. If subject's recon-all folder cannot be
 found, recon-all will be run. If not specified
 freesurfer data will be saved to {out\_dir}/freesurfer
 --models {Liem2016\_\_OCI\_norm,Liem2016\_\_full\_2samp\_training} [{Liem2016\_\_OCI\_norm,Liem2016\_\_full\_2samp\_training} ...]
 --skip\_missing Flag to skip not segmented subjects
 --license\_key LICENSE\_KEY
 FreeSurfer license key - letters and numbers after "\*"
 in the email you received after registration. To
 register (for free) visit
 https://surfer.nmr.mgh.harvard.edu/registration.html
 --n\_cpus N\_CPUS Number of CPUs/cores available to use.
 -v, --version show program's version number and exit

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. baracus.sh). For example:



```

#!/bin/bash
# sbatch --gres=lscratch:100 --mem=3g --cpus-per-task=48 --time=10:00:00 baracus.sh

set -e

function fail {
    echo "FAIL: $@" >&2
    exit 1  # signal failure
}

module load baracus
cp /data/$USER/BARACUS_TEST/ds001.tar.gz /lscratch/$SLURM_JOB_ID/.
cd /lscratch/$SLURM_JOB_ID
tar -xzvf ds001.tar.gz
baracus --participant_label 01 \
        --license_key /usr/local/apps/freesurfer/license.txt \
        --n_cpus $SLURM_CPUS_PER_TASK \
        /lscratch/$SLURM_JOB_ID/ds001 \
        /lscratch/$SLURM_JOB_ID \
        participant
cp -R /lscratch/$SLURM_JOB_ID/baracus /data/$USER/BARACUS_TEST/.
cp -R /lscratch/$SLURM_JOB_ID/freesurfer /data/$USER/BARACUS_TEST/.

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--gres=lscratch:#] [--cpus-per-task=#] [--mem=#] baracus.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. baracus.swarm). For example:



```

# swarm -f baracus.swarm --gres=lscratch:100 -g 4 -t 48 --time 10:00:00 --module baracus
cp /data/$USER/BARACUS_TEST/ds001.tar.gz /lscratch/$SLURM_JOB_ID/.; \
cd /lscratch/$SLURM_JOB_ID; \
tar -xzvf ds001.tar.gz; \
baracus --participant_label 01 \
        --license_key /usr/local/apps/freesurfer/license.txt \
        --n_cpus $SLURM_CPUS_PER_TASK \
        /lscratch/$SLURM_JOB_ID/ds001 \
        /lscratch/$SLURM_JOB_ID \
        participant; \
mv /lscratch/$SLURM_JOB_ID/baracus /data/$USER/BARACUS_TEST/baracus_01; \
mv /lscratch/$SLURM_JOB_ID/freesurfer /data/$USER/BARACUS_TEST/freesurfer_01
cp /data/$USER/BARACUS_TEST/ds001.tar.gz /lscratch/$SLURM_JOB_ID/.; \
cd /lscratch/$SLURM_JOB_ID; \
tar -xzvf ds001.tar.gz; \
baracus --participant_label 02 \
        --license_key /usr/local/apps/freesurfer/license.txt \
        --n_cpus $SLURM_CPUS_PER_TASK \
        /lscratch/$SLURM_JOB_ID/ds001 \
        /lscratch/$SLURM_JOB_ID \
        participant; \
mv /lscratch/$SLURM_JOB_ID/baracus /data/$USER/BARACUS_TEST/baracus_02; \
mv /lscratch/$SLURM_JOB_ID/freesurfer /data/$USER/BARACUS_TEST/freesurfer_02

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f baracus.swarm [--gres=lscratch:#] [-g #] [-t #] --module baracus
```

where


|  |  |  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- | --- | --- |
| -gres=lscratch:*#*  Number of Gigabytes of local disk space allocated per process (1 line in the swarm command file)
 | -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module baracus Loads the baracus module for each subjob in the swarm 
 | |
 | |
 | |
 | |








