

document.querySelector('title').textContent = 'MRIQC on Biowulf';
MRIQC on Biowulf


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



MRIQC (MRI quality control) is an application that automatically extracts image quality metrics from MRI scans and generates individual reports. MRIQC can fit a classifier to categorize datasets into "accept" or "exclude" categories.



### References:


* Esteban O, Birman D, Schaer M, Koyejo OO, Poldrack RA, Gorgolewski KJ.
 [*MRIQC: Advancing the automatic prediction of image quality in MRI from unseen sites.*](https://pubmed.ncbi.nlm.nih.gov/28945803/) 
 PLoS One. 2017 Sep 25;12(9):e0184661.


  

**IMPORTANT: (October 2021) The memory and cpu limit flags in mriqc do not work as intended.**   

The flags that are used to select an upper bound for mriqc memory processes and CPUs do not work as intended. If an mriqc job exceeds its memory allocation, the job will hang and D processes will be generated in the compute node where the job is running. We recommend carefully profiling memory consumption for a given type of dataset by running small-scale jobs, then gradually increasing the number of jobs when specific memory requirements have been established. Regarding the CPU flag not working, please use the `--exclusive` flag if using `sbatch` or `-t auto` if using swarm. Please get in touch with HPC staff if you need help with profiling memory requirements or allocating CPUs.



Documentation
* [MRIQC Documentation](https://mriqc.readthedocs.io)
* [Github page](https://github.com/poldracklab/mriqc)


Important Notes
* Module Name: mriqc (see [the modules page](/apps/modules.html) for more information)
* mriqc is a scientific pipeline with the potential to overload Biowulf's central filesystem. To avoid filesystem issues we recommend the following:
	1. Limit I/O with central Biowulf directories (e.g., /data, /home) by making use of local disk (/lscratch). You can make use of lscratch to store temporary working directories by directly assigning the temporary work directory of mriqc (-w option) to /lscratch/$SLURM\_JOB\_ID (remember to allocate enough space in /lscratch).
	* Estimate memory usage. Determine how much memory your jobs will need by benchmarking (see point 5 below). Something important to remember regarding memory allocation is to try to allocate some extra memory for your jobs as described in the resources described in point 5 below. 
	* Allocate all CPUs in a node. The flag that limits the number of CPUs in MRIQC do not work so please use the `--exclusive` flag if using `sbatch` or `-t auto` if using swarm.
	* Opt out of uploading mriqc metrics to mriqc website. You can disable the submission with the --no-sub flag.
	* Profile/benchmark mriqc jobs: We recommend making sure a given mriqc job can scale before launching a large number of jobs. You can do this by profiling/benchmarking small jobs (e.g., a swarm of 3 mriqc commands), then monitoring the jobs by using either the [user dashboard](https://hpc.nih.gov/dashboard/) or the commands jobhist, sjobs, squeue, and jobload (see [biowulf utilities](https://hpc.nih.gov/docs/biowulf_tools.html)). There are resources prepared by the HPC staff that go over how to monitor your jobs for benchmarking and profiling ([video](https://youtu.be/fLMJ8-t5bm4) and [slides](https://hpc.nih.gov/training/handouts/Effective_batch_system.pdf)). Once you have profiled a swarm with a few jobs and you have refined the memory, CPU and walltime requirements, it is probably safe to (gradually) increase the number of mriqc jobs. For many analysis pipelines one has no way of knowing in advance how much memory will be actually required in an HPC environment. This is why it is very important to profile/benchmark.



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

[user@cn3144 ~]$ **module load mriqc**

[user@cn3144 ~]$ **mriqc -h**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. mriqc.sh). For example (using example data):



```

#!/bin/bash
#SBATCH --job-name=mriqc
#SBATCH --gres=lscratch:20
#SBATCH --exclusive
#SBATCH --mem=16g
#SBATCH --time=72:00:00

module load mriqc/0.16.1

tar -C /lscratch/${SLURM_JOB_ID} -xf /usr/local/apps/mriqc/TEST_DATA/ds001.tar.gz

mriqc /lscratch/${SLURM_JOB_ID}/ds001 /lscratch/${SLURM_JOB_ID}/mriqc.out.ds001 \
      participant --participant_label sub-01 -w /lscratch/${SLURM_JOB_ID} \
      --no-sub 

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch mriqc.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. mriqc.swarm). For example:



```

mriqc /data/${USER}/BIDS-dataset/ds001/ /data/$USER/BIDS-dataset/mriqc.out.ds001 \
      participant --participant_label sub-01 -w /lscratch/${SLURM_JOB_ID} --no-sub \
mriqc /data/${USER}/BIDS-dataset/ds002/ /data/$USER/BIDS-dataset/mriqc.out.ds002 \
      participant --participant_label sub-02 -w /lscratch/${SLURM_JOB_ID} --no-sub \
mriqc /data/${USER}/BIDS-dataset/ds003/ /data/$USER/BIDS-dataset/mriqc.out.ds003 \
      participant --participant_label sub-03 -w /lscratch/${SLURM_JOB_ID} --no-sub \

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f mriqc.swarm [--gres=lscratch:#] [-g #] -t auto --module mriqc
```

where


|  |  |  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- | --- | --- |
| -gres=lscratch:*#*  Number of Gigabytes of local disk space allocated per process (1 line in the swarm command file)
 | -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file). We set this to `auto` to allocate all CPUs in each node.
 | --module mriqc Loads the mriqc module for each subjob in the swarm 
 | |
 | |
 | |
 | |








