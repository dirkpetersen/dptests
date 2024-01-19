

document.querySelector('title').textContent = 'Svtyper on HPC';
Svtyper on HPC


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


 Svtyper is a Bayesian genotyper for structural variants. 


Example dataset can be copied from */usr/local/apps/svtyper/testdata* 
 directory


Documentation
* <https://github.com/hall-lab/svtyper>



Important Notes
* Module Name: svtyper (see [the 
 modules page](/apps/modules.html) for more information)
* Example files in */usr/local/apps/svtyper/testdata*





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

[user@cn3144 ~]$ **module load svtyper**
[user@cn3144 ~]$ **svtyper -B NA12878.20.bam -S NA12878.20.splitters.bam > NA12878.20.gt.vcf**

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
module load svtyper
svtyper -B NA12878.20.bam -S NA12878.20.splitters.bam > NA12878.20.gt.vcf
```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--mem=#] batch.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. job.swarm). For example:



```

cd dir1; svtyper -B NA12878.20.bam -S NA12878.20.splitters.bam > NA12878.20.gt.vcf
cd dir2; svtyper -B NA12878.20.bam -S NA12878.20.splitters.bam > NA12878.20.gt.vcf
cd dir3; svtyper -B NA12878.20.bam -S NA12878.20.splitters.bam > NA12878.20.gt.vcf

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f job.swarm [-g #] --module svtyper
```

where
 

|  |  |
| --- | --- |
| -g *#*  | Number of Gigabytes of memory required for each process 
 (1 line in the swarm command file)  |
| --module  | Loads the module for each subjob in the swarm  |




