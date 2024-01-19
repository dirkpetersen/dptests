

document.querySelector('title').textContent = 'Rilseq on HPC';
Rilseq on HPC


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

  This package can be used to analyzed RILseq experiments. It is written 
 for a prokaryotic genome, without splice junction mapping and with some 
 additional features. The package handles the different stages processing 
 fastq files to pairs of interacting RNAs and some statistics. It does not 
 handle quality issues, adapter removing etc. so the fastq files should be 
 treated with cutadapt or equivalent before applying this package.


### References:

 * <http://www.cell.com/molecular-cell/fulltext/S1097-2765(16)30413-0>


Documentation * [https://github.com/asafpr/RILseq](https://github.com/asafpr/RILseq/)



Important Notes * Module Name: rilseq (see [the modules 
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

[user@cn3144 ~]$ **module load rilseq**
[user@cn3144 ~]$ **rilseq command**

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
module load rilseq
rilseq command
```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--mem=#] batch.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. job.swarm). For example:



```
cd dir1;rilseq command 1; rilseq command 2  
cd dir2;rilseq command 1; rilseq command 2  
cd dir3;rilseq command 1; rilseq command 2
```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f job.swarm [-t #] --module rilseq
```

where
 

|  |  |
| --- | --- |
| -g *#*  | Number of Gigabytes of memory required for each process 
 (1 line in the swarm command file)  |
| --module  | Loads the module for each subjob in the swarm  |




