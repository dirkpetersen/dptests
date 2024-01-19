

document.querySelector('title').textContent = 'Cnvnator on HPC';
Cnvnator on HPC


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


CNVnator is a tool for CNV discovery and genotyping from depth of read 
 mapping. 


 


Paper
* <https://www.ncbi.nlm.nih.gov/pubmed/21324876>



Documentation
* <https://github.com/abyzovlab/CNVnator>



Important Notes
* Module Name: cnvnator(see [the 
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

[user@cn3144 ~]$ **module load cnvnator**
[user@cn3144 ~]$ **cnvnator -root Output.root -chrom 1 2**
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
module load cnvnator
cnvnator -root Output.root -tree Input.bam
```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch batch.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. job.swarm). For example:



```

cd dir1; cnvnator -root Output.root -tree Input.bam
cd dir2; cnvnator -root Output.root -tree Input.bam
cd dir3; cnvnator -root Output.root -tree Input.bam

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f job.swarm --module cnvnator
```

where
 

|  |  |
| --- | --- |
| -g *#*  | Number of Gigabytes of memory required for each process 
 (1 line in the swarm command file)  |
| --module  | Loads the module for each subjob in the swarm  |




