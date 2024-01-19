

document.querySelector('title').textContent = 'Preseq on HPC';
Preseq on HPC


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

  The preseq package is aimed at predicting the yield of distinct reads from 
 a genomic library from an initial sequencing experiment. The estimates can 
 then be used to examine the utility of further sequencing, optimize the 
 sequencing depth, or to screen multiple libraries to avoid low complexity 
 samples.

 ### 


Documentation * <https://github.com/smithlabcode/preseq>



Important Notes * Module Name: preseq (see [the modules 
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

[user@cn3144 ~]$ **module load preseq**
[user@cn3144 ~]$ **preseq lc\_extrap -o yield\_estimates.txt input.bed**

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
module load preseq
preseq lc_extrap -o yield_estimates.txt input.bed
```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--mem=#] batch.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. job.swarm). For example:



```

cd dir1; preseq lc_extrap -o yield_estimates.txt input.bed
cd dir2; preseq lc_extrap -o yield_estimates.txt input.bed
cd dir3; preseq lc_extrap -o yield_estimates.txt input.bed

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f job.swarm [-g #]--module preseq
```

where
 

|  |  |
| --- | --- |
| -g *#*  | Number of Gigabytes of memory required for each process 
 (1 line in the swarm command file)  |
| --module  | Loads the module for each subjob in the swarm  |




