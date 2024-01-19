

document.querySelector('title').textContent = "gfatools";
gfatools on Biowulf


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



 gfatools is a set of tools for manipulating sequence graphs in the GFA or the [rGFA format](https://github.com/lh3/gfatools/blob/master/doc/rGFA.md).
 It has implemented parsing, subgraph and conversion to FASTA/BED.




Documentation
* [gfatools Main Site](https://github.com/lh3/gfatools)


Important Notes
* Module Name: gfatools (see [the modules page](/apps/modules.html) for more information)
* Environment variables set 
	+ GFATOOLS\_HOME* Example files in $GFATOOLS\_HOME/test



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

[user@cn3144 ~]$ **module load gfatools**

[user@cn3144 ~]$ **gfatools view -l MTh4502 -r 1 $GFATOOLS\_HOME/test/MT.gfa > sub.gfa**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. gfatools.sh). For example:



```

#!/bin/bash
set -e
module load gfatools
gfatools view -l MTh4502 -r 1 $GFATOOLS_HOME/test/MT.gfa > sub.gfa

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] gfatools.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. gfatools.swarm). For example:



```

gfatools stat sample1.gfa > sample1.txt
gfatools stat sample2.gfa > sample2.txt
gfatools stat sample3.gfa > sample3.txt
gfatools stat sample4.gfa > sample4.txt

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f gfatools.swarm [-g #] [-t #] --module gfatools
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module gfatools Loads the gfatools module for each subjob in the swarm 
 | |
 | |
 | |








