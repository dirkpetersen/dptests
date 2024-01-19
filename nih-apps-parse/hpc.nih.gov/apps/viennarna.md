

document.querySelector('title').textContent = 'Viennarna on HPC';
Viennarna on HPC


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

  The ViennaRNA Package consists of a C code library and several stand-alone 
 programs for the prediction and comparison of RNA secondary structures.


### References:

 * <http://dx.doi.org/10.1186/1748-7188-6-26>


Documentation * <https://www.tbi.univie.ac.at/RNA/#>



Important Notes * Module Name: viennarna (see [the 
 modules page](/apps/modules.html) for more information)





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

[user@cn3144 ~]$ **module load viennarna**
[user@cn3144 ~]$ **RNAfold -h**

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
module load viennarna
RNAfold < INPUT.fa

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch batch.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. job.swarm). For example:



```

cd dir1; RNAfold < INPUT.fa
cd dir2; RNAfold < INPUT.fa
cd dir3; RNAfold < INPUT.fa

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f job.swarm --module viennarna
```

where
 

|  |  |
| --- | --- |
| --module  | Loads the module for each subjob in the swarm  |




