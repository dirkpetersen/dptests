

document.querySelector('title').textContent = 'Control-FREEC on HPC';
Control-FREEC on HPC


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

  **Control-FREEC** is a tool for detection of copy-number 
 changes and allelic imbalances (including LOH) using deep-sequencing data 
 developed by the [Bioinformatics Laboratory 
 of Institut Curie (Paris)](http://bioinfo-out.curie.fr/).


Documentation * <https://github.com/BoevaLab/FREEC>



Important Notes
* Module Name: freec (see [the modules 
 page](/apps/modules.html) for more information)
* Example files in /usr/local/apps/freec/test.zip





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

[user@cn3144 ~]$ **module load freec**
[user@cn3144 ~]$ **freec -conf config\_ctrl.txt**

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
module load freec
freec -conf config_ctrl.txt
```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch batch.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. job.swarm). For example:



```

cd dir1; freec -conf config_ctrl.txt
cd dir2; freec -conf config_ctrl.txt
cd dir3; freec -conf config_ctrl.txt

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f job.swarm --module freec
```

where
 

|  |  |
| --- | --- |
| -g *#*  | Number of Gigabytes of memory required for each process 
 (1 line in the swarm command file)  |
| --module  | Loads the module for each subjob in the swarm  |




