

document.querySelector('title').textContent = 'NGMLR on Biowulf';
NGMLR on Biowulf


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



CoNvex Gap-cost alignMents for Long Reads (ngmlr) is a long-read mapper designed to sensitively align PacBilo or Oxford Nanopore to (large) reference genomes. It was designed to quickly and correctly align the reads, including those spanning (complex) structural variations. Ngmlr uses an SV aware k-mer search to find approximate mapping locations for a read and then a banded Smith-Waterman alignment algorithm to compute the final alignment. Ngmlr uses a convex gap cost model that penalizes gap extensions for longer gaps less than for shorter ones to compute precise alignments. The gap model allows ngmlr to account for both the sequencing error and real genomic variations at the same time and makes it especially effective at more precisely identifying the position of breakpoints stemming from structural variations. The k-mer search helps to detect and split reads that cannot be aligned linearly, enabling ngmlr to reliably align reads to a wide range of different structural variations including nested SVs (e.g. inversions flanked by deletions).



### References:


* Accurate detection of complex structural variations using single-molecule sequencing.
Fritz J. Sedlazeck, Philipp Rescheneder, Moritz Smolka, Han Fang, Maria Nattestad, Arndt von Haeseler & Michael C. Schatz.
Nature Methods
volume 15, 461–468 (2018).
[[html]](https://www.nature.com/articles/s41592-018-0001-7)


Documentation
* [NGMLR Main Site](https://github.com/philres/ngmlr)


Important Notes
* Module Name: ngmlr (see [the modules page](/apps/modules.html) for more information)
* Example files in /fdb/ngmlr



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive --cpus-per-task 4**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load ngmlr**

[user@cn3144 ~]$ **ngmlr -t $SLURM\_CPUS\_PER\_TASK -r reference.fasta -q reads.fastq -o test.sam** 

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. ngmlr.sh). For example:



```

#!/bin/sh
set -e
module load ngmlr

test -n "$SLURM_CPUS_PER_TASK" || SLURM_CPUS_PER_TASK=2

ngmlr -t $SLURM_CPUS_PER_TASK -r reference.fasta -q reads.fastq -o test.sam 

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] ngmlr.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. ngmlr.swarm). For example:



```

test -n "$SLURM_CPUS_PER_TASK" || SLURM_CPUS_PER_TASK=2; ngmlr -t $SLURM_CPUS_PER_TASK -r reference.fasta -q sample1.fastq -o sample1.sam
test -n "$SLURM_CPUS_PER_TASK" || SLURM_CPUS_PER_TASK=2; ngmlr -t $SLURM_CPUS_PER_TASK -r reference.fasta -q sample2.fastq -o sample2.sam
test -n "$SLURM_CPUS_PER_TASK" || SLURM_CPUS_PER_TASK=2; ngmlr -t $SLURM_CPUS_PER_TASK -r reference.fasta -q sample3.fastq -o sample3.sam
test -n "$SLURM_CPUS_PER_TASK" || SLURM_CPUS_PER_TASK=2; ngmlr -t $SLURM_CPUS_PER_TASK -r reference.fasta -q sample4.fastq -o sample4.sam

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f ngmlr.swarm [-g #] [-t #] --module ngmlr
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module ngmlr Loads the NGMLR module for each subjob in the swarm 
 | |
 | |
 | |








