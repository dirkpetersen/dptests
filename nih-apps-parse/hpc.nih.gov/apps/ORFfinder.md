

document.querySelector('title').textContent = 'ORFfinder on Biowulf';
ORFfinder on Biowulf


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



ORF finder searches for open reading frames (ORFs) in the DNA sequence you enter. The program returns the range of each ORF, along with its protein translation. Use ORF finder to search newly sequenced DNA for potential protein encoding segments, verify predicted protein using newly developed SMART BLAST or regular BLASTP.



Documentation
* [ORFfinder Main Site](https://www.ncbi.nlm.nih.gov/orffinder/)


Important Notes
* Module Name: ORFfinder (see [the modules page](/apps/modules.html) for more information)
* singlethreaded app



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

[user@cn3144 ~]$ module load ORFfinder

[user@cn3144 ~]$ ORFfinder -in /fdb/app_testdata/fasta/R64-1-1.cdna_nc.fa \
-id YPR161C -logfile text.out -out orf_test

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. ORFfinder.sh). For example:



```

#!/bin/bash
set -e
module load ORFfinder
ORFfinder -in /fdb/app_testdata/fasta/R64-1-1.cdna_nc.fa \
-id YPR161C -logfile text.out -out orf_test
```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch ORFfinder.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. ORFfinder.swarm). For example:



```

ORFfinder -in sample1.fa -out sample1 [...]
ORFfinder -in sample2.fa -out sample2 [...]
ORFfinder -in sample3.fa -out sample3 [...]
[...] rest of options

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f ORFfinder.swarm --module ORFfinder
```

where


|  |  |
| --- | --- |
| --module ORFfinder Loads the ORFfinder module for each subjob in the swarm 
 | |








