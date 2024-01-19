

document.querySelector('title').textContent = 'usearch on Biowulf';
usearch on Biowulf


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



USEARCH is a unique sequence analysis tool with thousands of users world-wide. USEARCH offers search and clustering algorithms that are often orders of magnitude faster than BLAST.



Documentation
* <http://www.drive5.com/usearch/>


Important Notes
* Module Name: usearch (see [the modules page](/apps/modules.html) for more information)
 * Single-threaded
 * Environment variables set:
	+ **usearch** -- path to usearch
	+ **USEARCH\_EXAMPLES** -- examples



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

[user@cn3144 ~]$ **module load usearch**
[user@cn3144 ~]$ **cp -R $USEARCH\_EXAMPLES/hmptut .**
[user@cn3144 ~]$ **cd hmptut/scripts**
[user@cn3144 ~]$ **bash run.bash**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. usearch.sh). For example:



```

#!/bin/bash
# --- This file is usearch.sh --
  
module load usearch
usearch -cluster_fast seqs.fasta -id 0.9 -centroids nr.fasta

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] usearch.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. usearch.swarm). For example:



```

usearch -cluster_fast seqs.fasta.1 -id 0.9 -centroids nr.fasta
usearch -cluster_fast seqs.fasta.2 -id 0.9 -centroids nr.fasta
usearch -cluster_fast seqs.fasta.3 -id 0.9 -centroids nr.fasta
usearch -cluster_fast seqs.fasta.4 -id 0.9 -centroids nr.fasta

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f usearch.swarm [-g #] [-t #] --module usearch
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module usearch Loads the usearch module for each subjob in the swarm 
 | |
 | |
 | |








