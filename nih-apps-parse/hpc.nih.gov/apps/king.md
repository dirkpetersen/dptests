

document.querySelector('title').textContent = 'KING on Biowulf';
KING on Biowulf


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



KING is a toolset to explore genotype data from a genome-wide association study (GWAS) or a sequencing project. KING can be used to check family relationship and flag pedigree errors by estimating kinship coefficients and inferring IBD segments for all pairwise relationships. Unrelated pairs can be precisely separated from close relatives with no false positives, with accuracy up to 3rd- or 4th-degree (depending on array or WGS) for --related and --ibdseg analyses, and up to 2nd-degree for --kinship analysis.



### References:


* Manichaikul A, Mychaleckyj JC, Rich SS, Daly K, Sale M, Chen WM.
 **[Robust relationship inference in genome-wide association studies.](https://www.ncbi.nlm.nih.gov/pubmed/20926424)**
Bioinformatics. 2010 Nov 15;26(22):2867-73.


Documentation
* [KING Main Site](http://people.virginia.edu/~wc9c/KING/)


Important Notes
* Module Name: king (see [the modules page](/apps/modules.html) for more information)
* Singlethreaded
* Environment variables set 
	+ KING\_EXAMPLES -- example files for KING



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

[user@cn3144 ~]$ **module load king**

[user@cn3144 ~]$ **cp $KING\_EXAMPLES/\* .**
[user@cn3144 ~]$ **king -b ex.bed --related**
[user@cn3144 ~]$ **king -b ex.bed --fam ex.fam --bim ex.bim --related**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. king.sh). For example:



```

#!/bin/bash
set -e
module load king
king -b my.bed --duplicate --cpus $SLURM_CPUS_PER_TASK
king -b my.bed --related --cpus $SLURM_CPUS_PER_TASK
king -b my.bed --related --degree --cpus $SLURM_CPUS_PER_TASK

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] king.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. king.swarm). For example:



```

king -b test1.bed --fam test1.fam --bim test1.bim --related
king -b test2.bed --fam test2.fam --bim test2.bim --related
king -b test3.bed --fam test3.fam --bim test3.bim --related
king -b test4.bed --fam test4.fam --bim test4.bim --related

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f king.swarm [-g #] [-t #] --module king
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module king Loads the king module for each subjob in the swarm 
 | |
 | |
 | |








