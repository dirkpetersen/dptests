

document.querySelector('title').textContent = 'Sicer on HPC';
Sicer on HPC


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


 Sicer is a clustering approach for identification of enriched domains from 
 histone modification ChIP-Seq data


### 


Documentation
* https://home.gwu.edu/~wpeng/Software.htm



Important Notes
* Module Name: SICER (see [the modules 
 page](/apps/modules.html) for more information)
* environment variables set
	+ SICERDIR
* Example files in /usr/local/apps/sicer/1.1/ex





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

[user@cn3144 ~]$ **module load sicer**
[user@cn3144 ~]$ **cd /data/$USER/sicer**
[user@cn3144 ~]$ **sh $SICERDIR/SICER.sh /data/$USER/sicer/ test.bed control.bed . hg18 1 200 150 0.74 600 .01**

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
module load sicer
sh $SICERDIR/SICER.sh /data/$USER/sicer/ test.bed control.bed . hg18 1 200 150 0.74 600 .01

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--mem=#] batch.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. job.swarm). For example:



```

cd dir1; sh $SICERDIR/SICER.sh /data/$USER/sicer/ test.bed control.bed . hg18 1 200 150 0.74 600 .01
cd dir2; sh $SICERDIR/SICER.sh /data/$USER/sicer/ test.bed control.bed . hg18 1 200 150 0.74 600 .01
cd dir3; sh $SICERDIR/SICER.sh /data/$USER/sicer/ test.bed control.bed . hg18 1 200 150 0.74 600 .01

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f job.swarm [-g #] --module sicer
```

where
 

|  |  |
| --- | --- |
| -g *#*  | Number of Gigabytes of memory required for each process 
 (1 line in the swarm command file)  |
| --module  | Loads the module for each subjob in the swarm  |




