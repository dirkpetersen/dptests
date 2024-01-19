

document.querySelector('title').textContent = 'VSEARCH on Biowulf';
VSEARCH on Biowulf


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



VSEARCH supports de novo and reference based chimera detection, clustering, full-length and prefix dereplication, rereplication, reverse complementation, masking, all-vs-all pairwise global alignment, exact and global alignment searching, shuffling, subsampling and sorting. It also supports FASTQ file analysis, filtering, conversion and merging of paired-end reads.



### References:


* Rognes T, Flouri T, Nichols B, Quince C, Mahé F.
 [**VSEARCH: a versatile open source tool for metagenomics.**](https://www.ncbi.nlm.nih.gov/pubmed/27781170)
*PeerJ. 2016 Oct 18;4*


Documentation
* [VSEARCH Main Site](https://github.com/torognes/vsearch)


Important Notes
* Module Name: vsearch (see [the modules page](/apps/modules.html) for more information)
* Multithreaded
* Example files in $VSEARCH\_EXAMPLES



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

[user@cn3144 ~]$ **module load vsearch**
[user@cn3144 ~]$ **ln -s /usr/local/apps/vsearch/vsearch-data/BioMarKs50k.fsa .**
[user@cn3144 ~]$ **vsearch --cluster\_fast BioMarKs50k.fsa --id 0.97 --centroids vsearch.out**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. vsearch.sh). For example:



```

#!/bin/bash
module load vsearch
vsearch --usearch_global queries.fsa --db database.fsa --id 0.9 --alnout alnout.txt

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] vsearch.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. vsearch.swarm). For example:



```

vsearch --usearch_global queries1.fsa --db database.fsa --id 0.9 --alnout alnout1.txt
vsearch --usearch_global queries2.fsa --db database.fsa --id 0.9 --alnout alnout2.txt
vsearch --usearch_global queries3.fsa --db database.fsa --id 0.9 --alnout alnout3.txt
vsearch --usearch_global queries4.fsa --db database.fsa --id 0.9 --alnout alnout4.txt

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f vsearch.swarm [-g #] [-t #] --module vsearch
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module vsearch Loads the VSEARCH module for each subjob in the swarm 
 | |
 | |
 | |








