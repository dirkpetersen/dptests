

document.querySelector('title').textContent = 'Samblaster on Biowulf';
Samblaster on Biowulf


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



samblaster is a program for marking duplicates and finding
discordant/split read pairs in read-id grouped paired-end SAM files. When
marking duplicates, samblaster will use about 20MB per 1M read pairs.

In a read-id grouped SAM file all alignments for a read-id (QNAME) are
continuous. Aligners naturally produce such files. They can also be created
by sorting a SAM file by read-id.



### References:


* Gregory G. Faust and Ira M. Hall, *SAMBLASTER: fast duplicate marking and
 structural variant read extraction*, Bioinformatics 2014, 30:
 2503-2505.
 [Pubmed](http://www.ncbi.nlm.nih.gov/pubmed/24812344) | 
 [PMC](http://www.ncbi.nlm.nih.gov/pmc/articles/pmid/24812344/) | 
 [Journal](http://bioinformatics.oxfordjournals.org/content/30/17/2503)


Documentation
* [Samblaster Main Site](https://github.com/GregoryFaust/samblaster)


Important Notes
* Module Name: samblaster (see [the modules page](/apps/modules.html) for more information)
* Singlethreaded app
* Example files in /usr/local/apps/samblaster/TEST\_DATA



Run samblaster on a bam file sorted by read name with duplicates already
marked. Save discordant pairs to `disc.sam` and split reads
to `split.sam`


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

[user@cn3144 ~]$ **module load samblaster**
[+] Loading samblaster 0.1.22

[user@cn3144 ~]$ **samtools view -h /usr/local/apps/samblaster/TEST\_DATA/test.bam \**
  **| samblaster --ignoreUnmated -a -e -d disc.sam -s split.sam -o /dev/null**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. samblaster.sh). For example:



```

#!/bin/bash

module load samtools samblaster || exit 1
samtools view -h /path/to/input.bam \
  | samblaster -e -d disc.sam -s split.sam -o /dev/null


```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch samblaster.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. samblaster.swarm). For example:



```

samtools view -h /path/to/input1.bam \
  | samblaster -e -d disc1.sam -s split1.sam -o /dev/null
samtools view -h /path/to/input2.bam \
  | samblaster -e -d disc2.sam -s split2.sam -o /dev/null
samtools view -h /path/to/input3.bam \
  | samblaster -e -d disc3.sam -s split3.sam -o /dev/null

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f samblaster.swarm --module samblaster
```

where


|  |  |
| --- | --- |
| --module samblaster Loads the samblaster module for each subjob in the swarm 
 | |








