

document.querySelector('title').textContent = 'BRASS on Biowulf';
BRASS on Biowulf


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



BRASS analyses one or more related BAM files of paired-end sequencing to determine potential rearrangement breakpoints.

There are several stages the main component being:
* Collect read-pairs where both ends map but NOT marked as properly-paired.
 * Perform grouping based on mapped locations
* Filter
* Run assembly
* Annotate with GRASS





Documentation
* [BRASS on GitHub](https://github.com/cancerit/BRASS)


Important Notes
* Module Name: BRASS (see [the modules page](/apps/modules.html) for more information)
* BRASS is part of the Cancer Genome Project and is closely related to the programs [Ascat NGS](/apps/ascatNgs.html) and [cgpBattenberg](/apps/cgpBattenberg.html) as well as the utilites [VAGrENT](https://github.com/cancerit/VAGrENT) and [PCAP-core](https://github.com/cancerit/PCAP-core). All of these programs can be added to your path using the cancerit-wgs module. To get the most recent versions of all of these, use the cancerit-wgs/latest module version.
* Multithreaded. Note: the BRASS documentation recommends a max of 2 CPUs during the 'input' phase.



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive --cpus-per-task=4 --mem=20g**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load BRASS**

[user@cn3144 ~]$  **brass.pl -c $SLURM\_CPUS\_PER\_TASK -o myout -t tumour.bam -n normal.bam**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. BRASS.sh). For example:



```

#!/bin/bash
set -e
module load BRASS
brass.pl -c $SLURM_CPUS_PER_TASK -o myout -t tumour.bam -n normal.bam

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=4] [--mem=20g] BRASS.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. BRASS.swarm). For example:



```

brass.pl -c $SLURM_CPUS_PER_TASK -o myout1 -t tumour1.bam -n normal.bam
brass.pl -c $SLURM_CPUS_PER_TASK -o myout2 -t tumour2.bam -n normal.bam
[...]

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f BRASS.swarm -g 15 -t 4 --module BRASS
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module BRASS Loads the BRASS module for each subjob in the swarm 
 | |
 | |
 | |








