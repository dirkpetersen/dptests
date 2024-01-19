

document.querySelector('title').textContent = 'Tractseg on Biowulf';
Tractseg on Biowulf


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






### References:


* Wasserthal J, Neher P, Maier-Hein KH
 [*TractSeg - Fast and accurate white matter tract segmentation*](https://pubmed.ncbi.nlm.nih.gov/30086412/) 
 Neuroimage. 2018 Dec;183:239-253.


  

Documentation
* [Github page](https://github.com/MIC-DKFZ/TractSeg)


Important Notes
* Module Name: tractseg (see [the modules page](/apps/modules.html) for more information)
* Environment variables: TRACTSEG\_HOME (application directory), TRACTSEG\_TESTDATA (test dataset).
* Limit I/O with central Biowulf directories (e.g., /data, /home) by making use of local disk (/lscratch).
* Estimate memory usage. Determine how much memory your jobs will need by benchmarking. Allocate some extra memory for your jobs as described in the resources described below. 
* Profile/benchmark tractseg jobs: We recommend making sure a given tractseg job can scale before launching a large number of jobs. You can do this by profiling/benchmarking small jobs (e.g., a swarm of 3 tractseg commands), then monitoring the jobs by using either the [user dashboard](https://hpcnihapps.cit.nih.gov/auth/dashboard/) or the commands jobhist, sjobs, squeue, and jobload (see [biowulf utilities](https://hpc.nih.gov/docs/biowulf_tools.html)). There are resources prepared by the HPC staff that go over how to monitor your jobs for benchmarking and profiling ([video](https://youtu.be/fLMJ8-t5bm4), and [slides](https://hpc.nih.gov/training/handouts/Effective_batch_system.pdf)). Once you have profiled a swarm with a few jobs and you have refined the memory, and walltime requirements, it is probably safe to (gradually) increase the number of tractseg jobs. For many analyses pipelines one has no way of knowing in advance how much memory will be actually required in an HPC environment. This is why it is very important to profile/benchmark.



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

[user@cn3144 ~]$ **module load tractseg**

[user@cn3144 ~]$ **TractSeg --help**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. tractseg.sh). For example (using test data below):



```

#!/bin/bash
#SBATCH --job-name=tractseg
#SBATCH --gres=lscratch:20
#SBATCH --exclusive
#SBATCH --mem=32g
#SBATCH --time=72:00:00

module load fmriprep/20.2.1

tar -C /lscratch/${SLURM_JOB_ID} -xf /usr/local/apps/fmriprep/TEST_DATA/ds001.tar.gz

fmriprep /lscratch/${SLURM_JOB_ID}/ds001 /lscratch/${SLURM_JOB_ID}/fmriprep.out.ds001 \
         participant --participant_label sub-01 -w /lscratch/${SLURM_JOB_ID} \
         --notrack --use-aroma

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch fmriprep.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. fmriprep.swarm). For example:



```

fmriprep /data/${USER}/BIDS-dataset/ds001/ /data/$USER/BIDS-dataset/fmriprep.out.ds001 \
         participant --participant_label sub-01 -w /lscratch/${SLURM_JOB_ID} \
         --use-aroma --notrack
fmriprep /data/${USER}/BIDS-dataset/ds002/ /data/$USER/BIDS-dataset/fmriprep.out.ds002 \
         participant --participant_label sub-02 -w /lscratch/${SLURM_JOB_ID} \
         --use-aroma --notrack
fmriprep /data/${USER}/BIDS-dataset/ds003/ /data/$USER/BIDS-dataset/fmriprep.out.ds003 \
         participant --participant_label sub-03 -w /lscratch/${SLURM_JOB_ID} \
         --use-aroma --notrack 

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f fmriprep.swarm [--gres=lscratch:#] [-g #] -t auto --module fmriprep
```

where


|  |  |  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- | --- | --- |
| -gres=lscratch:*#*  Number of Gigabytes of local disk space allocated per process (1 line in the swarm command file)
 | -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file). We set this to `auto` to allocate all CPUs in each node.
 | --module fmriprep Loads the fmriprep module for each subjob in the swarm 
 | |
 | |
 | |
 | |








