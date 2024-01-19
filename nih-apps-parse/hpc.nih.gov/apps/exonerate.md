

document.querySelector('title').textContent = 'Exonerate on Biowulf';
Exonerate on Biowulf


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



Exonerate is a generic tool for pairwise sequence comparison. It allows you to align sequences using a many alignment models, either exhaustive dynamic programming or a variety of heuristics.



### References:


* Slater GS and Birney E (2005) Automated generation of heuristics for biological sequence comparison. BMC Bioinformatics 6:31; doi: [10.1186/1471-2105-6-31](https://doi.org/10.1186/1471-2105-6-31)


Documentation
* [Exonerate Main Site](https://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate)


Important Notes
* Module Name: exonerate (see [the modules page](/apps/modules.html) for more information)



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

[user@cn3144 ~]$ **module load exonerate**
[+] Loading exonerate, version 2.2.0...
[user@cn3144 ~]$ **exonerate query.fasta target.fasta**
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. exonerate.sh). For example:



```

#!/bin/bash
set -e
module load exonerate
exonerate query.fasta target.fasta

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] exonerate.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. exonerate.swarm). For example:



```

exonerate query01.fasta target.fasta
exonerate query02.fasta target.fasta
exonerate query03.fasta target.fasta
exonerate query04.fasta target.fasta

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f exonerate.swarm [-g #] [-t #] --module exonerate
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module exonerate Loads the exonerate module for each subjob in the swarm 
 | |
 | |
 | |








