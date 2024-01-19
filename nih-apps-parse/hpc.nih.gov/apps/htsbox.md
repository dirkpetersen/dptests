

document.querySelector('title').textContent = "htsbox";
htsbox on Biowulf


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



HTSbox is a fork of early HTSlib. It is a collection of small experimental tools manipulating HTS-related files. While some of these tools are already part of the official SAMtools package, others are for niche use cases.




Documentation
* [htsbox Main Site](https://github.com/lh3/htsbox)


Important Notes
* Module Name: htsbox (see [the modules page](/apps/modules.html) for more information)
 * Multithreaded with the -t flag for certain subcommands.
* Environment variables set 
	+ HTSBOX\_HOME* Example files in $HTSBOX\_HOME/test



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

[user@cn3144 ~]$ **module load htsbox**

[user@cn3144 ~]$ **htsbox abreak -l 0 $HTSBOX\_HOME/test/ex3.sam**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. htsbox.sh). For example:



```

#!/bin/bash
set -e
module load htsbox
htsbox abreak -l 0 $HTSBOX_HOME/test/ex3.sam

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] htsbox.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. htsbox.swarm). For example:



```

htsbox qualbin -t $SLURM_CPUS_PER_TASK -bm7 in1.bam
htsbox qualbin -t $SLURM_CPUS_PER_TASK -bm7 in2.bam
htsbox qualbin -t $SLURM_CPUS_PER_TASK -bm7 in3.bam
htsbox qualbin -t $SLURM_CPUS_PER_TASK -bm7 in4.bam

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f htsbox.swarm [-g #] [-t #] --module htsbox
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module htsbox Loads the htsbox module for each subjob in the swarm
 | |
 | |
 | |








