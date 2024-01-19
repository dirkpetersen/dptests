

document.querySelector('title').textContent = 'bam2fastq on HPC';
bam2fastq on HPC


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

  BAM-format files are used to store alignment information and unaligned reads 
 from next-generation sequencing machines. This tool is intended to extract 
 raw sequences (with qualities) from a BAM file. 

 ### 


Documentation * <https://gsl.hudsonalpha.org/information/software/bam2fastq>


 
Important Notes


 * Module Name: bam2fastq (see [the 
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

[user@cn3144 ~]$ **module load bam2fastq**
[user@cn3144 ~]$ **bam2fastq input.bam -o outfile#**
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit). Create a batch input file (e.g. bam2fastq.sh). For example:



```

#!/bin/bash
set -e
module load bam2fastq
bam2fastq input.bam -o outfile#
```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch bam2fastq.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. bam2fastq.swarm). For example:



```

cd dir1; bam2fastq input.bam -o outfile#
cd dir2; bam2fastq input.bam -o outfile#
cd dir3; bam2fastq input.bam -o outfile#

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f bam2fastq.swarm [-t #] --module bam2fastq
```

where
 

|  |  |  |  |
| --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process 
 (1 line in the swarm command file) 
 | --module  Loads the module for each subjob in the swarm 
  | |
 | |








