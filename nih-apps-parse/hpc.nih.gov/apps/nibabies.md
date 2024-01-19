

document.querySelector('title').textContent = 'Nibabies on Biowulf';
Nibabies on Biowulf


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



Nibabies is an MRI preprocessing pipeline for the infant and neonate brain.



### References:


* Mathias Goncalves, Christopher Markiewicz, Martin Styner, Lucille Moore, Kathy Snider, Eric Earl, Christopher Smyser, 
 Lilla Zollei, Russell Poldrack, Oscar Esteban, Eric Feczko, Damien Fair.
 *NiBabies: A robust preprocessing workflow tailored for neonate and infant MRI.* 
 Organization for Human Brain Mapping, 2021 Annual Meeting. Abstract # 2452.


### Web site


* [Github page](https://github.com/nipreps/nibabies)


  

Documentation
* [Nibabies Documentation](https://https://github.com/nipreps/nibabies)


Important Notes
* nibabies on Biowulf is a singularity container built directly from a docker container. However, users do not need to execute singularity directly or bind directories manually because we have provided a wrapper shell script that takes care of that. 
* Module Name: nibabies (see [the modules page](/apps/modules.html) for more information)
* nibabies is a scientific pipeline with the potential to overload Biowulf's central filesystem. To avoid filesystem issues we recommend the following:
	1. Limit I/O with central Biowulf directories (e.g., /data, /home) by making use of local disk (/lscratch). You can make use of lscratch to store temporary working directories by directly assigning the temporary work directory of fmriprep (-w option) to /lscratch/$SLURM\_JOB\_ID (remember to allocate enough space in /lscratch).
	* Limit memory usage. Once you determine how much memory your jobs will need (by benchmarking, but see point 5 below), It is still good idea to limit the memory used by nibabies by using the option --mem\_mb. Something important to remember regarding memory allocation is to try to allocate some extra memory for your jobs as described in the resources described in point 6 below. 
	* Limit multi-threading. Biowulf's nibabies has been installed as a container and some of the applications in it will hyper-thread unless you explicitly limit the number of threads. You can limit the number of threads that nibabies is allowed to use across all processes by using the --nthreads flag.
	* Opt out of uploading nibabies metrics to nibabies website. You can disable the submission with the --notrack flag.
	* Make use of the flag --stop-on-first-crash flag to force a stop if issues occur.
	* Profile/benchmark nibabies jobs: We recommend making sure a given nibabies job can scale before launching a large number of jobs. You can do this by profiling/benchmarking small jobs (e.g., a swarm of 3 nibabies commands), then monitoring the jobs by using either the [user dashboard](https://hpcnihapps.cit.nih.gov/auth/dashboard/) or the commands jobhist, sjobs, squeue, and jobload (see [biowulf utilities](https://hpc.nih.gov/docs/biowulf_tools.html)). There are resources prepared by the HPC staff that go over how to monitor your jobs for benchmarking and profiling ([video](https://youtu.be/fLMJ8-t5bm4), and [slides](https://hpc.nih.gov/training/handouts/Effective_batch_system.pdf)). Once you have profiled a swarm with a few jobs and you have refined the memory, CPU (and walltime) requirements, it is probably safe to (gradually) increase the number of nibabies jobs. For many analyses pipelines one has no way of knowing in advance how much memory or CPUs will be actually required in an HPC environment. This is why it is very important to profile/benchmark.



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int), run the program and download some imaging datasets to test. Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load nibabies**

[user@cn3144 ~]$ **nibabies --help**

[user@cn3144 ~]$ **cd /data/$USER**

[user@cn3144 ~]$ **mkdir NIBABIES\_TEST**

[user@cn3144 ~]$ **cd NIBABIES\_TEST**

[user@cn3144 ~]$ **module load aws**

[user@cn3144 ~]$ **aws s3 sync --no-sign-request s3://openneuro.org/ds003778 ds003778**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. nibabies.sh). For example:



```

#!/bin/bash
#SBATCH --job-name=nibabies
#SBATCH --gres=lscratch:20
#SBATCH --exclusive
#SBATCH --mem=32g
#SBATCH --time=72:00:00

module load nibabies

nibabies /data/$USER/NIBABIES_TEST/ds003778 /data/$USER/NIBABIES_TEST/nibabies.out.ds003778 \
         participant --participant-label sub-01 -w /lscratch/${SLURM_JOB_ID} \
         --age-months 12 --notrack

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch nibabies.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. nibabies.swarm). For example:



```

nibabies /data/$USER/NIBABIES_TEST/ds003778 /data/$USER/NIBABIES_TEST/nibabies.out.ds003778-01 \
         participant --participant-label sub-01 -w /lscratch/${SLURM_JOB_ID} \
         --age-months 12 --notrack
nibabies /data/$USER/NIBABIES_TEST/ds003778 /data/$USER/NIBABIES_TEST/nibabies.out.ds003778-02 \
         participant --participant-label sub-02 -w /lscratch/${SLURM_JOB_ID} \
         --age-months 12 --notrack
nibabies /data/$USER/NIBABIES_TEST/ds003778 /data/$USER/NIBABIES_TEST/nibabies.out.ds003778-03 \
         participant --participant-label sub-03 -w /lscratch/${SLURM_JOB_ID} \
         --age-months 12 --notrack

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f nibabies.swarm [--gres=lscratch:#] [-g #] [-t #] --module nibabies
```

where


|  |  |  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- | --- | --- |
| -gres=lscratch:*#*  Number of Gigabytes of local disk space allocated per process (1 line in the swarm command file)
 | -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module nibabies Loads the nibabies module for each subjob in the swarm 
 | |
 | |
 | |
 | |








