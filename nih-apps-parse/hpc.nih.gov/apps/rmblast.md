

document.querySelector('title').textContent = 'RMBlast on Biowulf';
RMBlast on Biowulf


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



RMBlast is a [RepeatMasker](repeatmasker.html)-compatible version of the standard NCBI blastn program.
RMBlast supports RepeatMasker searches by adding a few necessary features to the stock NCBI blastn program. These include:

* Support for custom matrices ( without KA-Statistics ).
* Support for cross\_match-like complexity adjusted scoring. Cross\_match is Phil Green's seeded smith-waterman search algorithm.
* Support for cross\_match-like masklevel filtering.






Documentation
* [RMBlast Main Site](http://www.repeatmasker.org/RMBlast.html)


Important Notes
* Module Name: rmblast (see [the modules page](/apps/modules.html) for more information)



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **rmblastn**
[user@cn3144 ~]$ 
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. rmblast.sh). For example:



```

#!/bin/sh
set -e
module load rmblast
rmblastn -query query.fasta -subject reference.fasta

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] rmblast.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. rmblast.swarm). For example:



```

rmblastn -query query01.fasta -subject reference.fasta
rmblastn -query query02.fasta -subject reference.fasta
rmblastn -query query03.fasta -subject reference.fasta
rmblastn -query query04.fasta -subject reference.fasta

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f rmblast.swarm [-g #] [-t #] --module rmblast
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module rmblast Loads the rmblast module for each subjob in the swarm 
 | |
 | |
 | |








