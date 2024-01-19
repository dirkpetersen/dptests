

document.querySelector('title').textContent = "gutSMASH";
gutSMASH on Biowulf


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



 gutSMASH is a tool that has been developed to systematically evaluate the metabolic potential of anaerobic bacteria in the gut by predicting both known and novel anaerobic metabolic gene clusters (MGCs) from the gut microbiome.
 The gutSMASH detection rules have been validated using a curated dataset.



### References:


* Andreu, V.P., Augustijn, H.E., Chen, L., Zhernakova, A., Fu, J., Fischbach, M.A., Dodd, D. and Medema, M.H.
 [**A systematic analysis of metabolic pathways in the human gut microbiota.**](https://doi.org/10.1101/2021.02.25.432841)
*bioRxiv, pp.2021-02..*


Documentation
* [gutSMASH Main Site](https://gutsmash.bioinformatics.nl/help.html)
* [gutSMASH on GitHub](https://github.com/victoriapascal/gutsmash)


Important Notes
This application requires a [graphical connection using NX](/docs/connect.html#nx)


* Module Name: gutsmash (see [the modules page](/apps/modules.html) for more information)
 * Multithreaded: set the -c/--cpus flag to match your allocation.
 * This application produces HTML reports. You can use [hpcdrive](/docs/hpcdrive.html) to view these reports on your local workstation.
* Environment variables set 
	+ GUTSMASH\_HOME* Example files in $GUTSMASH\_HOME/test* Reference data in /fdb/gutsmash/



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive --cpus-per-task 4**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load gutsmash**

[user@cn3144 ~]$ **run\_gutsmash.py --cpus $SLURM\_CPUS\_PER\_TASK --minimal $GUTSMASH\_HOME/test/Streptomyces\_coelicolor.gbk**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. gutsmash.sh). For example:



```

#!/bin/bash
set -e
module load gutsmash
run_gutsmash --cpus $SLURM_CPUS_PER_TASK --minimal $GUTSMASH_HOME/test/Streptomyces_coelicolor.gbk

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] gutsmash.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. gutsmash.swarm). For example:



```

run_gutsmash.py --cpus $SLURM_CPUS_PER_TASK --minimal sample1.gbk
run_gutsmash.py --cpus $SLURM_CPUS_PER_TASK --minimal sample2.gbk
run_gutsmash.py --cpus $SLURM_CPUS_PER_TASK --minimal sample3.gbk

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f gutsmash.swarm [-g #] [-t #] --module gutsmash
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module gutsmash Loads the gutSMASH module for each subjob in the swarm 
 | |
 | |
 | |








