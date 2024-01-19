

document.querySelector('title').textContent = 'breseq on Biowulf';
breseq on Biowulf


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




From the breseq documentation:




> 
> 
>  breseq is a computational pipeline for the analysis of short-read
>  re-sequencing data (e.g. Illumina, 454, IonTorrent, etc.). It uses
>  reference-based alignment approaches to predict mutations in a sample
>  relative to an already sequenced genome. breseq is intended for microbial
>  genomes (<10 Mb) and re-sequenced samples that are only slightly diverged
>  from the reference sequence (<1 mutation per 1000 bp).
> 
>  breseq‘s primary advantages over other software programs are that it can:
>  1. Accurately predict new sequence junctions, such as those associated with 
>  mobile element insertions.
> 2. Integrate multiple sources of evidence for genetic changes into mutation predictions.
> 3. Produce annotated output describing biologically relevant mutational events.
> 
> 
> 
>  breseq was initially developed to analyze data from the Lenski long-term evolution 
>  experiment with E. coli. However, breseq may be generally useful to researchers who are:
> 
>  1. Tracking mutations over time in microbial evolution experiments.
> 2. Checking strains for unwanted second-site mutations after genetic manipulations.
> 3. Identifying mutations that occur during strain improvement or after 
>  long-term culture of engineered strains.
> 4. Discovering what mutations arise in pathogens during infection or cause antibiotic resistance.
> 
> 
> 


### References:


* D. E. Deatherage and J. E. Barrick. 
 *Identification of mutations in laboratory-evolved microbes from next-generation sequencing 
 data using breseq*. Methods Mol. Biol. 1151: 165–188 (2014)
 [PubMed](https://www.ncbi.nlm.nih.gov/pubmed/24838886) | 
 [PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4239701/) | 
 [Journal](https://link.springer.com/protocol/10.1007%2F978-1-4939-0554-6_12)



Documentation
* breseq [manual](http://barricklab.org/twiki/pub/Lab/ToolsBacterialGenomeResequencing/documentation/index.html)
* breseq [GitHub repo](https://github.com/barricklab/breseq)


Important Notes
* Module Name: breseq (see [the modules page](/apps/modules.html) for more information)
* breseq can use multiple threads. Please allocate one CPU for each thread
* Example files in `BRESEQ_TEST_DATA`



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. 
In this sample session we will analyze an E. coli strain that evolved for 20000 generations
in the [long term evolution experiment](http://myxo.css.msu.edu/ecoli/).



```

[user@biowulf]$ **sinteractive --mem=5g --cpus-per-task=4 --gres=lscratch:10**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144]$ **cd /lscratch/$SLURM\_JOB\_ID**
[user@cn3144]$ **module load breseq**
[user@cn3144]$ **cp -L ${BRESEQ\_TEST\_DATA:-none}/\* .**
[user@cn3144]$ **ls -lh**
total 275M
-rw-r--r-- 1 user group  11M Aug  8 13:46 NC_012967.gbk
-rw-r--r-- 1 user group 133M Aug  8 13:46 SRR030257_1.fastq.gz
-rw-r--r-- 1 user group 132M Aug  8 13:46 SRR030257_2.fastq.gz

[user@cn3144]$ **breseq -j $SLURM\_CPUS\_PER\_TASK -r NC\_012967.gbk \
 SRR030257\_1.fastq.gz SRR030257\_2.fastq.gz**
...
---> bowtie2  :: version 2.3.4.1 [/usr/local/apps/bowtie/2-2.3.4.1/bin/bowtie2]
---> R        :: version 3.5.0 [/usr/local/apps/R/3.5/3.5.0_build2/bin/R]
+++   NOW PROCESSING Read and reference sequence file input
  READ FILE::SRR030257_1
...
[user@cn3144]$ **cp -r output data /path/to/where/you/would/like/the/output**
[user@cn3144]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf]$

```

The output directory contains summary html files.



Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. breseq.sh), which uses the input file 'breseq.in'. For example:



```

#!/bin/bash
module load breseq/0.36.1 || exit 1
wd=$PWD
cd /lscratch/$SLURM_JOB_ID || exit 1
cp -L $BRESEQ_TEST_DATA/* .
breseq -j $SLURM_CPUS_PER_TASK -r NC_012967.gbk \
    SRR030257_1.fastq.gz SRR030257_2.fastq.gz
cp -r output $wd

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --gres=lscratch:10 --cpus-per-task=6 --mem=5g breseq.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. breseq.swarm). For example:



```

breseq -j $SLURM_CPUS_PER_TASK -r ref.gbk sample1_reads.fastq.gz
breseq -j $SLURM_CPUS_PER_TASK -r ref.gbk sample2_reads.fastq.gz
breseq -j $SLURM_CPUS_PER_TASK -r ref.gbk sample3_reads_R1.fastq.gz sample3_reads_R2.fastq.gz

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f breseq.swarm -g 5 -t 4 --module breseq/0.33.0
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g #  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t #  Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module breseq  Loads the breseq module for each subjob in the swarm 
 | |
 | |
 | |








