

document.querySelector('title').textContent = 'Bamtools on Biowulf';
Bamtools on Biowulf


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

  BamTools provides both a programmer's API and an end-user's toolkit for handling 
 BAM files. It was developed by  [Derek Barnett in the Marth lab at Boston College](https://github.com/pezmaster31/bamtools/wiki).


### References:

 * <https://academic.oup.com/bioinformatics/article/27/12/1691/255399>


Documentation * <https://github.com/pezmaster31/bamtools/wiki>

 * Module Name: bamtools (see [the 
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

[user@cn3144 ~]$ **module load bamtools**
[user@cn3144 ~]$ **bamtools convert -format fastq -in in.bam -out out.fastq**
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit). Create a batch input file (e.g. bamtools.sh). For example:



```

#!/bin/bash
set -e
module load bamtools
bamtools convert -format fastq -in in.bam -out out.fastq
```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--mem=#] bamtools.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. bamtools.swarm). For example:



```

cd dir1; bamtools convert -format fastq -in in.bam -out out.fastq
cd dir2; bamtools convert -format fastq -in in.bam -out out.fastq
cd dir3; bamtools convert -format fastq -in in.bam -out out.fastq

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f bamtools.swarm [-g #] --module bamtools
```

where
 

|  |  |  |  |
| --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in 
 the swarm command file) 
 | --module bamtools Loads the bamtools module for each subjob in the swarm 
  | |
 | |








