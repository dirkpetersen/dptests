

document.querySelector('title').textContent = 'Smriprep on Biowulf';
Smriprep on Biowulf


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



Smriprep is a computational pipeline for preprocessing of structural MRI datasets. The output of smriprep can be used as an input to [dmriprep](https://hpc.nih.gov/apps/dmriprep.html) and [fmriprep](https://hpc.nih.gov/apps/fmriprep.html).



Documentation
* [Smriprep Documentation](https://www.nipreps.org/smriprep/)
* [Github page](https://github.com/nipreps/smriprep)


Important Notes
* Module Name: smriprep (see [the modules page](/apps/modules.html) for more information)
* smriprep is a scientific pipeline with the potential to overload Biowulf's central filesystem. To avoid filesystem issues we recommend the following:
	1. Limit I/O with central Biowulf directories (e.g., /data, /home) by making use of local disk (/lscratch). You can make use of lscratch to store temporary working directories by directly assigning the temporary work directory of smriprep (-w option) to /lscratch/$SLURM\_JOB\_ID (remember to allocate enough space in /lscratch).
	* Estimate memory usage. Determine how much memory your jobs will need by benchmarking (see point 3 below). Allocate some extra memory for your jobs as described in the resources described in point 3 below. 
	* Profile/benchmark smriprep jobs: We recommend making sure a given smriprep job can scale before launching a large number of jobs. You can do this by profiling/benchmarking small jobs (e.g., a swarm of 3 smriprep commands), then monitoring the jobs by using either the [user dashboard](https://hpc.nih.gov/dashboard/) or the commands jobhist, sjobs, squeue, and jobload (see [biowulf utilities](https://hpc.nih.gov/docs/biowulf_tools.html)). There are resources prepared by the HPC staff that go over how to monitor your jobs for benchmarking and profiling ([video](https://youtu.be/fLMJ8-t5bm4), and [slides](https://hpc.nih.gov/training/handouts/Effective_batch_system.pdf)). Once you have profiled a swarm with a few jobs and you have refined the memory, and walltime requirements, it is probably safe to (gradually) increase the number of smriprep jobs. For many analyses pipelines one has no way of knowing in advance how much memory will be actually required in an HPC environment. This is why it is very important to profile/benchmark.
	* Opt out of uploading fmriprep metrics to fmriprep website. You can disable the submission with the --notrack flag.



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

[user@cn3144 ~]$ **module load smriprep**

[user@cn3144 ~]$ **smriprep -h**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. smriprep.sh). For example (using test data below):



```

#!/bin/bash
#SBATCH --job-name=smriprep
#SBATCH --gres=lscratch:20
#SBATCH --exclusive
#SBATCH --mem=32g
#SBATCH --time=72:00:00

module load smriprep

smriprep /lscratch/${SLURM_JOB_ID}/input001 /lscratch/${SLURM_JOB_ID}/smriprep.out.ds001 \
         participant --participant_label sub-01 -w /lscratch/${SLURM_JOB_ID}

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch smriprep.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. smriprep.swarm). For example:



```

smriprep /data/${USER}/BIDS-dataset/input001/ /data/$USER/BIDS-dataset/smriprep.out.ds001 \
         participant --participant_label sub-01 -w /lscratch/${SLURM_JOB_ID} \
smriprep /data/${USER}/BIDS-dataset/input002/ /data/$USER/BIDS-dataset/smriprep.out.ds002 \
         participant --participant_label sub-02 -w /lscratch/${SLURM_JOB_ID} \
smriprep /data/${USER}/BIDS-dataset/input003/ /data/$USER/BIDS-dataset/smriprep.out.ds003 \
         participant --participant_label sub-03 -w /lscratch/${SLURM_JOB_ID} \

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f smriprep.swarm [--gres=lscratch:#] [-g #] -t auto --module smriprep
```

where


|  |  |  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- | --- | --- |
| -gres=lscratch:*#*  Number of Gigabytes of local disk space allocated per process (1 line in the swarm command file)
 | -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file). We set this to `auto` to allocate all CPUs in each node.
 | --module smriprep Loads the smriprep module for each subjob in the swarm 
 | |
 | |
 | |
 | |








