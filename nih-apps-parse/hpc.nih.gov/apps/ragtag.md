

document.querySelector('title').textContent = 'RagTag on Biowulf';
RagTag on Biowulf


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


RagTag is a collection of software tools for scaffolding and improving modern genome assemblies.


Tasks include:


* Homology-based misassembly correction
* Homology-based assembly scaffolding and patching
* Scaffold merging


Ragtag also provides command line utilities for working with common genome assembly file formats.


### References:


* Alonge, Michael, et al. 
 [**RaGOO: fast and accurate reference-guided scaffolding of draft genomes.**](https://pubmed.ncbi.nlm.nih.gov/31661016/)
*Genome biology 20.1 (2019): 1-17.*


Documentation
* [RagTag Main Site](https://github.com/malonge/RagTag/wiki)


Important Notes
* Module Name: ragtag (see [the modules page](/apps/modules.html) for more information)
 * singlethreaded
 * Reference data in /data/genome/fasta/



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

[user@cn3144 ~]$ **module load ragtag**

[user@cn3144 ~]$ **ragtag.py --help**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. ragtag.sh). For example:



```

#!/bin/bash
set -e
module load ragtag
ragtag.py correct ref.fasta query.fasta

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] ragtag.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. ragtag.swarm). For example:



```

ragtag.py correct ref.fasta query1.fasta
ragtag.py correct ref.fasta query2.fasta
ragtag.py correct ref.fasta query3.fasta
ragtag.py correct ref.fasta query4.fasta

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f ragtag.swarm [-g #] [-t #] --module ragtag
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module ragtag Loads the ragtag module for each subjob in the swarm 
 | |
 | |
 | |








