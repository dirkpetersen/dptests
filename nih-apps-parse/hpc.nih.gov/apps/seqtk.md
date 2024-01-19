

document.querySelector('title').textContent = 'seqtk on Biowulf';
seqtk on Biowulf


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

  Seqtk is a fast and lightweight tool for processing sequences in the FASTA 
 or FASTQ format. It seamlessly parses both FASTA and FASTQ files which can 
 also be optionally compressed by gzip. 

 ### 


Documentation * <https://github.com/lh3/seqtk>


Important Notes * Module Name: seqtk (see [the modules 
 page](/apps/modules.html) for more information)



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

[user@cn3144 ~]$ **module load seqtk**[user@cn3144 ~]$ **seqtk seq -a in.fq.gz > out.fa**
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. seqtk.sh). For example:



```

#!/bin/bash
set -e
module load seqtk
seqtk seq -a in.fq.gz > out.fa
```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--mem=#] seqtk.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. seqtk.swarm). For example:



```

cd dir1;seqtk seq -a in.fq.gz > out.fa
cd dir2;seqtk seq -a in.fq.gz > out.fa
cd dir3;seqtk seq -a in.fq.gz > out.fa

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f seqtk.swarm [-g #] --module seqtk
```

where
 

|  |  |  |  |
| --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in 
 the swarm command file) 
 | --module seqtk Loads the seqtk module for each subjob in the swarm 
  | |
 | |








