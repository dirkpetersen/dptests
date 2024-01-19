

document.querySelector('title').textContent = 'PRSice on Biowulf';
PRSice on Biowulf


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



PRSice (pronounced 'precise') is a Polygenic Risk Score software for calculating, applying, evaluating and plotting the results of polygenic risk scores (PRS) analyses.



### References:


* PRSice: Polygenic Risk Score software, Euesden, Lewis, O'Reilly, Bioinformatics (2015) 31 (9):1466-1468. doi: [10.1093/bioinformatics/btu848](https://doi.org/10.1093/bioinformatics/btu848)


Documentation
* [PRSice Main Site](https://choishingwan.github.io/PRSice/)
* [User Support Group](https://groups.google.com/forum/#!forum/prsice)


Important Notes
* Module Name: prsice (see [the modules page](/apps/modules.html) for more information)
* Environment variables set 
	+ PRSICE\_HOME



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
[user@cn3144 ~]$ cd /data/$USER
[user@cn3144 ~]$ mkdir TEST
[user@cn3144 ~]$ cd TEST
[user@cn3144 ~] cp /usr/local/apps/prsice/TEST_DATA/* .
[user@cn3144 ~]$ **module load prsice**
[user@cn3144 ~]$ **PRSice.R \
 --prsice $(which PRSice) \
 --base TOY\_BASE\_GWAS.assoc \
 --target TOY\_TARGET\_DATA \
 --thread 1 \
 --stat OR \
 --binary-target T**
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. prsice.sh). For example:



```

#!/bin/sh
set -e
module load prsice

test -n "$SLURM_CPUS_PER_TASK" || SLURM_CPUS_PER_TASK=1
PRSice.R \
 --prsice $(which PRSice) \
 --base TOY_BASE_GWAS.assoc \
 --target TOY_TARGET_DATA \
 --thread $SLURM_CPUS_PER_TASK \
 --stat OR
 --binary-target T

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] prsice.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. prsice.swarm). For example:



```

PRSice.R \
 --prsice $(which PRSice) \
 --base TOY_BASE_GWAS.assoc \
 --target TARGET_DATA_1 \
 --thread $SLURM_CPUS_PER_TASK \
 --stat OR
 --binary-target T
PRSice.R \
 --prsice $(which PRSice) \
 --base TOY_BASE_GWAS.assoc \
 --target TARGET_DATA_2 \
 --thread $SLURM_CPUS_PER_TASK \
 --stat OR
 --binary-target T
PRSice.R \
 --prsice $(which PRSice) \
 --base TOY_BASE_GWAS.assoc \
 --target TARGET_DATA_3 \
 --thread $SLURM_CPUS_PER_TASK \
 --stat OR
 --binary-target T

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f prsice.swarm [-g #] -t # --module prsice
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module prsice Loads the PRSice module for each subjob in the swarm 
 | |
 | |
 | |








