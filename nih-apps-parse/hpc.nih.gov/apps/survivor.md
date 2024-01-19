

document.querySelector('title').textContent = 'SURVIVOR on Biowulf';
SURVIVOR on Biowulf


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



SURVIVOR is a tool set for simulating/evaluating SVs, merging and comparing SVs within and among samples, and includes various methods to reformat or summarize SVs.



### References:


* Transient structural variations have strong effects on quantitative traits and reproductive isolation in fission yeast.
Jeffares, Daniel C; Jolly, Clemency; Hoti, Mimoza; Speed, Doug; Shaw, Liam; Rallis, Charalampos; Balloux, Francois; Dessimoz, Christophe; Bähler, Jürg; Sedlazeck, Fritz J.
Nature communications, Vol. 8, 14061, 24.01.2017, p. 1-11. DOI:[10.1038/NCOMMS14061](https://doi.org/10.1038/NCOMMS14061)


Documentation
* [SURVIVOR Main Site](https://github.com/fritzsedlazeck/SURVIVOR)
* [SURVIVOR wiki](https://github.com/fritzsedlazeck/SURVIVOR/wiki)


Important Notes
* Module Name: survivor (see [the modules page](/apps/modules.html) for more information)



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

[user@cn3144 ~]$ **module load survivor**

[user@cn3144 ~]$ **ls \*vcf > sample-files**
[user@cn3144 ~]$ **SURVIVOR merge sample-files 1000 2 1 1 0 30 sample-merged.vcf**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. survivor.sh). For example:



```

#!/bin/sh
set -e
module load survivor
SURVIVOR bincov lowMQ.cov 10 2 > lowMQ.bed

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] survivor.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. survivor.swarm). For example:



```

SURVIVOR merge sample1-files 1000 2 1 1 0 30 sample1-merged.vcf
SURVIVOR merge sample2-files 1000 2 1 1 0 30 sample2-merged.vcf
SURVIVOR merge sample3-files 1000 2 1 1 0 30 sample3-merged.vcf
SURVIVOR merge sample4-files 1000 2 1 1 0 30 sample4-merged.vcf

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f survivor.swarm [-g #] [-t #] --module survivor
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module survivor Loads the SURVIVOR module for each subjob in the swarm 
 | |
 | |
 | |








