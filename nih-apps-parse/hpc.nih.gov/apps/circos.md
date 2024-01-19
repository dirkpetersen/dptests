

document.querySelector('title').textContent = 'Circos on HPC';
Circos on HPC


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


 Circos is a program for the generation of publication-quality, circularly 
 composited renditions of genomic data and related annotations. Circos is 
 particularly suited for visualizing alignments, conservation and intra and 
 inter-chromosomal relationships. Also, Circos is useful to visualize any 
 type of information that benefits from a circular layout. Thus, although 
 it has been designed for the field of genomics, it is sufficiently flexible 
 to be used in other data domains. 


### References:


* [Krzywinski, 
 M. et al. Circos: an Information Aesthetic for Comparative Genomics. Genome 
 Res (2009) 19:1639-1645](http://genome.cshlp.org/content/early/2009/06/15/gr.092759.109.abstract)


Documentation
* <http://circos.ca/tutorials/lessons/>



Important Notes * Module Name: circos (see [the modules 
 page](/apps/modules.html) for more information)





Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load circos**
[user@cn3144 ~]$ **circos commands...**
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```




Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. batch.sh). For example:



```

#!/bin/bash
set -e
module load circos
circos commands
```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--mem=#] batch.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. job.swarm). For example:



```

cd dir1; circos commnads...
cd dir2; circos commnads...
cd dir3; circos commnads...
cd dir4; circos commnads...

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f job.swarm [-g #] --module circos
```

where
 

|  |  |
| --- | --- |
| -g *#*  | Number of Gigabytes of memory required for each process 
 (1 line in the swarm command file)  |
| --module  | Loads the module for each subjob in the swarm  |




