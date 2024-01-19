

document.querySelector('title').textContent = "LAST";
LAST on Biowulf


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



 LAST is designed for moderately large data (e.g. genomes, DNA reads,
 proteomes). It's especially geared toward:
 * Finding rearrangements and recombinations (last-split)
* Finding DNA-versus-protein related regions, especially protein fossils.
* Unusual data, e.g. AT-rich DNA, because it can fit parameters to the data and calculate significance.
* Sensitive DNA-DNA search, due to fitting, sensitive seeding, and
 calculating significance.


 It can also indicate the confidence/uncertainty of each column in an
 alignment, and use sequence quality data in a rigorous fashion.



### References:


* Kiełbasa SM, Wan R, Sato K, Horton P, Frith MC.
 [**Adaptive seeds tame genomic sequence comparison.**](http://genome.cshlp.org/content/21/3/487.long)
*Genome Res. 2011 21(3):487-93.*
* See [**Detailed papers on LAST**](https://gitlab.com/mcfrith/last/-/blob/main/doc/last-papers.rst) for a full list.


Documentation
* [LAST Main Site](https://gitlab.com/mcfrith/last)
* [LAST Cookbook](https://gitlab.com/mcfrith/last/-/blob/main/doc/last-cookbook.rst)
* [Documentation page listing](https://gitlab.com/mcfrith/last/-/tree/main/doc)


Important Notes
* Module Name: last (see [the modules page](/apps/modules.html) for more information)
 * Multithreaded using the -P flag on lastal* Environment variables set 
	+ LAST\_HOME* Example files in $LAST\_HOME/examples



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

[user@cn3144 ~]$ **module load last**

[user@cn3144 ~]$ **lastdb humdb $LAST\_HOME/examples/humanMito.fa**
[user@cn3144 ~]$ **lastal humdb $LAST\_HOME/examples/fuguMito.fa > myalns.maf**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. last.sh). For example:



```

#!/bin/bash
set -e
module load last
lastdb humdb $LAST_HOME/examples/humanMito.fa
lastal humdb $LAST_HOME/examples/fuguMito.fa > myalns.maf

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] last.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. last.swarm). For example:



```

lastal -P $SLURM_CPUS_PER_TASK humdb sample1.fa > s1.maf
lastal -P $SLURM_CPUS_PER_TASK humdb sample2.fa > s2.maf
lastal -P $SLURM_CPUS_PER_TASK humdb sample3.fa > s3.maf
lastal -P $SLURM_CPUS_PER_TASK humdb sample4.fa > s4.maf

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f last.swarm [-g #] [-t #] --module last
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module last Loads the last module for each subjob in the swarm 
 | |
 | |
 | |








