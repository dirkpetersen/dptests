

document.querySelector('title').textContent = 'modkit on HPC';
modkit on HPC


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


 bioinformatics tool for working with modified bases from Oxford Nanopore. Specifically for converting modBAM to bedMethyl files using best practices, but also manipulating modBAM files and generating summary statistics.



### 


Documentation
* <https://nanoporetech.github.io/modkit/>


 
Important Notes



* Module Name: modkit (see [the 
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

[user@cn3144 ~]$ **module load modkit**
[user@cn3144 ~]$ **modkit pileup input\_reads.bam output/path/pileup.bed --log-filepath pileup.log**
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit). Create a batch input file (e.g. modkit.sh). For example:



```

#!/bin/bash
set -e
module load modkit
modkit pileup input_reads.bam output/path/pileup.bed --log-filepath pileup.log

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
[user@biowulf ~]$ **sbatch modkit.sh**
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. modkit.swarm). For example:



```

modkit pileup dir1/input_reads.bam dir1/pileup.bed --log-filepath pileup1.log
modkit pileup dir2/input_reads.bam dir2/pileup.bed --log-filepath pileup2.log
modkit pileup dir3/input_reads.bam dir3/pileup.bed --log-filepath pileup3.log

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f modkit.swarm [-g #] [-t #] --module modkit
```

where
 

|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process
 (1 line in the swarm command file)
 || -t *#*  Number of threads allocated for each process
 (1 line in the swarm command file)
 | --module  Loads the module for each subjob in the swarm
  | |
 | |

 | |








