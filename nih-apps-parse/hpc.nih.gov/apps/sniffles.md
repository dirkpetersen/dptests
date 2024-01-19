

document.querySelector('title').textContent = 'sniffles on Biowulf';
sniffles on Biowulf


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



Sniffles is a structural variation caller using third generation sequencing (PacBio or Oxford Nanopore). It detects all types of SVs (10bp+) using evidence from split-read alignments, high-mismatch regions, and coverage analysis. Please note the current version of Sniffles requires sorted output from [BWA-MEM](bwa.html) (use -M and -x parameter) or [NGMLR](ngmlr.html) with the optional SAM attributes enabled



### References:


* Accurate detection of complex structural variations using single-molecule sequencing.
Fritz J. Sedlazeck, Philipp Rescheneder, Moritz Smolka, Han Fang, Maria Nattestad, Arndt von Haeseler & Michael C. Schatz.
Nature Methods
volume 15, 461–468 (2018).
[[html]](https://www.nature.com/articles/s41592-018-0001-7)


Documentation
* [sniffles Main Site](https://github.com/fritzsedlazeck/Sniffles)


Important Notes
* Module Name: sniffles (see [the modules page](/apps/modules.html) for more information)
* Example files in /fdb/sniffles* sniffles uses 4 threads by default.
Specify --threads and allocate the corresponding number of threads with your job submission.



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

[user@cn3144 ~]$ **module load sniffles**
[+] Loading sniffles, version 2.0.2...
[user@cn3144 ~]$ **sniffles -i /fdb/sniffles/reads\_region.bam -v test.vcf**
Estimating parameter...
Max dist between aln events: 4
Max diff in window: 50
Min score ratio: 2
Avg DEL ratio: 0.0432235
Avg INS ratio: 0.0740733
Start parsing... 21
Switch Chr 21
Finalizing  ..
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. sniffles.sh). For example:



```

#!/bin/sh
set -e
module load sniffles
sniffles -i /fdb/sniffles/reads_region.bam -v test.vcf

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] sniffles.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. sniffles.swarm). For example:



```

sniffles -i sample1.bam -v sample1.vcf --threads $SLURM_CPUS_PER_TASK
sniffles -i sample2.bam -v sample2.vcf --threads $SLURM_CPUS_PER_TASK
sniffles -i sample3.bam -v sample3.vcf --threads $SLURM_CPUS_PER_TASK
sniffles -i sample4.bam -v sample4.vcf --threads $SLURM_CPUS_PER_TASK

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f sniffles.swarm [-g #] -t 4 --module sniffles
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module sniffles Loads the sniffles module for each subjob in the swarm 
 | |
 | |
 | |








