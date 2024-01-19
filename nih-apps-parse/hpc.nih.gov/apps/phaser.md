

document.querySelector('title').textContent = 'phASER on Biowulf';
phASER on Biowulf


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



phASER stands for *phasing and Allele Specific Expression from RNA-seq*. It performs haplotype phasing using read alignments in BAM format from both DNA and RNA based assays, and provides measures of haplotypic expression for RNA based assays.



### References:


* [Castel, Stephane E., et al. "Rare variant phasing and haplotypic expression from RNA sequencing with phASER." *Nature communications* 7.1 (2016): 1-6.](https://www.nature.com/articles/ncomms12817)


Documentation
* [phASER on GitHub](https://github.com/secastel/phaser)
* [phASER tutorial](https://stephanecastel.wordpress.com/2017/02/15/how-to-generate-ase-data-with-phaser/)


Important Notes
* Module Name: phaser (see [the modules page](/apps/modules.html) for more information)
 * phASER has a --threads option. Testing shows that --cpus-per-task should be set to twice the number specified by --threads.
 * Example data can be found in /usr/local/apps/phaser/1.1.1/testdata* Do **not** preface script calls with python. Just call the phaser\*.sh scripts directly and they will use the correct version of Python.



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf ~]$ **sinteractive -c8 --mem=10g --gres=lscratch:20**
salloc.exe: Pending job allocation 63023261
salloc.exe: job 63023261 queued and waiting for resources
salloc.exe: job 63023261 has been allocated resources
salloc.exe: Granted job allocation 63023261
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn0873 are ready for job
srun: error: x11: no local DISPLAY defined, skipping

[user@cn0873 ~]$ **cd /lscratch/$SLURM\_JOB\_ID**

[user@cn0873 63023261]$ **cp -r /usr/local/apps/phaser/1.1.1/testdata .**

[user@cn0873 63023261]$ **cd testdata/**

[user@cn0873 testdata]$ **module load phaser**
[+] Loading phaser  1.1.1  on cn0873
[+] Loading singularity  3.6.1  on cn0873

[user@cn0873 testdata]$ **phaser.py --vcf NA06986.vcf.gz \
 --bam NA06986.2.M\_111215\_4.bam --paired\_end 1 --mapq 255 --baseq 10 \
 --sample NA06986 --blacklist hg19\_hla.bed \
 --haplo\_count\_blacklist hg19\_haplo\_count\_blacklist.bed --threads 4 \
 --o phaser\_test\_case**

##################################################
              Welcome to phASER v1.1.1
  Author: Stephane Castel (scastel@nygenome.org)
  Updated by: Bishwa K. Giri (bkgiri@uncg.edu)
##################################################

Completed the check of dependencies and input files availability...

STARTED "Read backed phasing and ASE/haplotype analyses" ...
    DATE, TIME : 2020-08-14, 14:14:10
[...snip]
     COMPLETED using 1176416 reads in 481 seconds using 4 threads
     PHASED  23919 of 2142443 all variants (= 0.011164) with at least one other variant
     GENOME WIDE PHASE CORRECTED  1 of 2142443 variants (= 0.000000)
     Global maximum memory usage: 2822.19 (mb)

COMPLETED "Read backed phasing" of sample NA06986 in 00:08:31 hh:mm:ss
DATE, TIME : 2020-08-14, 14:22:42

The End.

[user@cn0873 testdata]$ **phaser\_gene\_ae.py \
 --haplotypic\_counts phaser\_test\_case.haplotypic\_counts.txt \
 --features gencode.v19.GRCh37.genes.bed --o phaser\_test\_case\_gene\_ae.txt**

##################################################
          Welcome to phASER Gene AE v1.2.0
  Author: Stephane Castel (scastel@nygenome.org)
##################################################

#1 Loading features...
#2 Loading haplotype counts...
sys:1: DtypeWarning: Columns (5,16,17) have mixed types. Specify dtype option on import or set low_memory=False.
#3 Processing results...
    BAM: NA06986.2.M_111215_4
          generating feature level haplotypic counts...
          outputting feature haplotype counts...

[user@cn0873 testdata]$ **exit**
exit
srun: error: cn0873: task 0: Exited with exit code 130
salloc.exe: Relinquishing job allocation 63023261
salloc.exe: Job allocation 63023261 has been revoked.

[user@biowulf ~]$


```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. phaser.sh). For example:



```

#!/bin/bash
set -e
module load phaser
phaser.py --vcf NA06986.vcf.gz --bam NA06986.2.M_111215_4.bam --paired_end 1 \
    --mapq 255 --baseq 10  --sample NA06986 --blacklist hg19_hla.bed \
    --haplo_count_blacklist hg19_haplo_count_blacklist.bed --threads 4 \
    --o phaser_test_case

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] phaser.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. phaser.swarm). For example:



```

phaser_gene_ae.py \
    --haplotypic_counts phaser_test_case.haplotypic_counts1.txt \
    --features gencode.v19.GRCh37.genes.bed --o phaser_test_case_gene_ae1.txt
phaser_gene_ae.py \
    --haplotypic_counts phaser_test_case.haplotypic_counts2.txt \
    --features gencode.v19.GRCh37.genes.bed --o phaser_test_case_gene_ae2.txt
phaser_gene_ae.py \
    --haplotypic_counts phaser_test_case.haplotypic_counts3.txt \
    --features gencode.v19.GRCh37.genes.bed --o phaser_test_case_gene_ae3.txt
phaser_gene_ae.py \
    --haplotypic_counts phaser_test_case.haplotypic_counts4.txt \
    --features gencode.v19.GRCh37.genes.bed --o phaser_test_case_gene_ae4.txt

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f phaser.swarm [-g #] [-t #] --module phaser
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module phaser Loads the phaser module for each subjob in the swarm 
 | |
 | |
 | |








