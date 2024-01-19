

document.querySelector('title').textContent = 'EDD on Biowulf';
EDD on Biowulf


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



EDD is a ChIP-seq peak caller for detection of megabase domains of enrichment.



### References:


* E Lund, AR Oldenburg, and P Collas
 [Enriched domain detector: a program for detection of wide genomic enrichment domains robust against local variations](https://academic.oup.com/nar/article/42/11/e92/1432049)
*Nucl. Acids Res. 42 (2014)*


Documentation
* [EDD on GitHub](https://github.com/CollasLab/edd)


Important Notes
* Module Name: edd (see [the modules page](/apps/modules.html) for more information)
 * Multithreaded
 * Environment variables set:
	+ **EDD\_EXAMPLES**



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive --cpus-per-task=8**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load edd**
[user@cn3144 ~]$ **ln -s $EDD\_EXAMPLES/hg19.sizes .**
[user@cn3144 ~]$ **ln -s $EDD\_EXAMPLES/\*.bam .**
[user@cn3144 ~]$ **touch unalign.bed**
[user@cn3144 ~]$ **edd -p $SLURM\_CPUS\_PER\_TASK -n 10000 hg19.sizes unalign.bed SRX447385.bam SRX447386.bam output\_dir**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. edd.sh). For example:



```

#!/bin/bash

# Generate hg19 chrom sizes
module load ucsc
fetchChromSizes hg19 > hg19.sizes
module unload ucsc

# Create empty file for unalignable regions
touch unalign_empty

# Run the command, collecing output in output directory
ml edd
edd -p $SLURM_CPUS_PER_TASK -n 10000 hg19.sizes unalign.bed chip_seq.bam input.bam output_dir

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] edd.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. edd.swarm). For example:



```

edd -p $SLURM_CPUS_PER_TASK hg19.sizes unalign.bed chip_seq_1.bam input_1.bam output_dir_1
edd -p $SLURM_CPUS_PER_TASK hg19.sizes unalign.bed chip_seq_2.bam input_2.bam output_dir_2
edd -p $SLURM_CPUS_PER_TASK hg19.sizes unalign.bed chip_seq_3.bam input_3.bam output_dir_3
edd -p $SLURM_CPUS_PER_TASK hg19.sizes unalign.bed chip_seq_4.bam input_4.bam output_dir_4

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f edd.swarm [-g #] [-t #] --module edd
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module edd Loads the edd module for each subjob in the swarm 
 | |
 | |
 | |








