

document.querySelector('title').textContent = "tandem-genotypes";
tandem-genotypes on Biowulf


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



tandem-genotypes finds changes in length of tandem repeats, from "long" DNA reads aligned to a genome.



### References:


* Mitsuhashi, S., Frith, M.C., Mizuguchi, T. et al.
 [**Tandem-genotypes: robust detection of tandem repeat expansions from long DNA reads.**](https://doi.org/10.1186/s13059-019-1667-6)
*Genome Biol 20, 58 (2019).*


Documentation
* [tandem-genotypes Main Site](https://github.com/mcfrith/tandem-genotypes)


Important Notes
* Module Name: tandem-genotypes (see [the modules page](/apps/modules.html) for more information)
* This application produces plots via the tandem-genotypes-plot command. You can use [hpcdrive](/docs/hpcdrive.html) to view them on your local workstation.
* Environment variables set 
	+ TANDEMGENOTYPES\_HOME* Example files in $TANDEMGENOTYPES\_HOME/tests



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

[user@cn3144 ~]$ **module load tandem-genotypes**

[user@cn3144 ~]$ **tandem-genotypes $TANDEMGENOTYPES\_HOME/tests/microsat.txt $TANDEMGENOTYPES\_HOME/tests/nano.maf**
# tandem-genotypes /usr/local/apps/tandem-genotypes/1.9.0/tests/microsat.txt /usr/local/apps/tandem-genotypes/1.9.0/tests/nano.maf
chr22	41994883	41994923	TG	.	.	-7,-6,-4,-3,-1,0,0,2,3,4,8	-16,-14,-12,-11,-11,-11,-9,-8,-7
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. tandem-genotypes.sh). For example:



```

#!/bin/bash
set -e
module load tandem-genotypes
tandem-genotypes $TANDEMGENOTYPES_HOME/tests/microsat.txt $TANDEMGENOTYPES_HOME/tests/nano.maf > tg.txt

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] tandem-genotypes.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. tandem-genotypes.swarm). For example:



```

tandem-genotypes $TANDEMGENOTYPES_HOME/tests/microsat.txt $TANDEMGENOTYPES_HOME/tests/sample1.maf > tg1.txt
tandem-genotypes $TANDEMGENOTYPES_HOME/tests/microsat.txt $TANDEMGENOTYPES_HOME/tests/sample2.maf > tg2.txt
tandem-genotypes $TANDEMGENOTYPES_HOME/tests/microsat.txt $TANDEMGENOTYPES_HOME/tests/sample3.maf > tg3.txt
tandem-genotypes $TANDEMGENOTYPES_HOME/tests/microsat.txt $TANDEMGENOTYPES_HOME/tests/sample4.maf > tg4.txt

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f tandem-genotypes.swarm [-g #] [-t #] --module tandem-genotypes
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module tandem-genotypes Loads the tandem-genotypes module for each subjob in the swarm 
 | |
 | |
 | |








