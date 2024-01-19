

document.querySelector('title').textContent = 'Mixcr on Biowulf';
Mixcr on Biowulf


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



MiXCR is a universal software for fast and accurate analysis of T- and B- cell receptor repertoire sequencing data.
Description

Documentation
* <https://mixcr.readthedocs.org/en/latest/>


Important Notes
* Module Name: mixcr (see [the modules page](/apps/modules.html) for more information)
 * Multithreaded



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive --cpus-per-task=4 --mem=10g**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load mixcr**

[user@cn3144 ~]$ **mixcr align -r log.txt -t $SLURM\_CPUS\_PER\_TASK \
 /data/$USER/S1.fastq.gz \
 /data/$USER/S2.fastq.gz \
 aln.vdjca** 

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. mixcr.sh). For example:



```

#!/bin/bash
set -e
module load mixcr
mixcr align -r log.txt -t $SLURM_CPUS_PER_TASK \
		/data/$USER/S1.fastq.gz \
		/data/$USER/S2.fastq.gz \
		aln.vdjca

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=4 --mem=10g mixcr.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. mixcr.swarm). For example:



```

cd /data/$USER/dir1; mixcr align -r log.txt -t $SLURM_CPUS_PER_TASK S1.fastq.gz S2.fastq.gz aln.vdjca
cd /data/$USER/dir2; mixcr align -r log.txt -t $SLURM_CPUS_PER_TASK S1.fastq.gz S2.fastq.gz aln.vdjca
cd /data/$USER/dir3; mixcr align -r log.txt -t $SLURM_CPUS_PER_TASK S1.fastq.gz S2.fastq.gz aln.vdjca
cd /data/$USER/dir4; mixcr align -r log.txt -t $SLURM_CPUS_PER_TASK S1.fastq.gz S2.fastq.gz aln.vdjca

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f mixcr.swarm -g 10 -t 4 --module mixcr
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module mixcr Loads the mixcr module for each subjob in the swarm 
 | |
 | |
 | |








