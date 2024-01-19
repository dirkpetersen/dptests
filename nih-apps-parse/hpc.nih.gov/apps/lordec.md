

document.querySelector('title').textContent = 'LoRDEC on Biowulf';
LoRDEC on Biowulf


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



LoRDEC processes data coming from high throughput sequencing machines of the second and third generations. These data are called sequencing reads, or simply reads for short. Technically speaking it processes short reads and long reads to correct errors in the long reads.



### References:


* Leena Salmela & Eric Rivals 
 [**LoRDEC: accurate and efficient long read error correction**](https://pubmed.ncbi.nlm.nih.gov/25165095/)
*Bioinformatics (2014), Dec 15;30(24):3506-14*


Documentation
* [LoRDEC Main Site](https://www.lirmm.fr/~rivals/lordec/)


Important Notes
This application requires a [graphical connection using NX](/docs/connect.html#nx)


* Module Name: lordec (see [the modules page](/apps/modules.html) for more information)
* Some commands are multithreaded



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

[user@cn3144 ~]$ **module load lordec**

[user@cn3144 ~]$ **lordec-correct -2 illumina.fasta -k 19 -s 3 -i pacbio.fasta -o pacbio-corrected.fasta**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. lordec.sh). For example:



```

#!/bin/bash
set -e
module load lordec
lordec-trim -i corrected.reads.fasta -o trimmed.reads.fasta

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] lordec.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. lordec.swarm). For example:



```

lordec-correct -2 illumina.fasta -k 19 -s 3 -i pacbio_1.fasta -o pacbio_1-corrected.fasta
lordec-correct -2 illumina.fasta -k 19 -s 3 -i pacbio_2.fasta -o pacbio_2-corrected.fasta
lordec-correct -2 illumina.fasta -k 19 -s 3 -i pacbio_3.fasta -o pacbio_3-corrected.fasta
lordec-correct -2 illumina.fasta -k 19 -s 3 -i pacbio_4.fasta -o pacbio_4-corrected.fasta

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f lordec.swarm [-g #] [-t #] --module lordec
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module lordec Loads the lordec module for each subjob in the swarm 
 | |
 | |
 | |








