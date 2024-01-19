

document.querySelector('title').textContent = 'QTLtools on Biowulf';
QTLtools on Biowulf


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



QTLtools is a tool set for molecular QTL discovery and analysis. It allows to go from the raw sequence data to collection of molecular Quantitative Trait Loci (QTLs) in few easy-to-perform steps.



### Reference:


* [Delaneau, Olivier, et al. "A complete tool set for molecular QTL discovery and analysis." *Nature communications* 8.1 (2017): 1-7.](https://www.nature.com/articles/ncomms15452)


Documentation
* [QTLtools documentation](https://qtltools.github.io/qtltools/)
* [QTLtools on GitHub](https://github.com/qtltools/qtltools)


Important Notes
* Module Name: qtltools (see [the modules page](/apps/modules.html) for more information)



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf ~]$ **sinteractive -c4 --mem=8g --gres=lscratch:10**
salloc.exe: Pending job allocation 61635379
salloc.exe: job 61635379 queued and waiting for resources
salloc.exe: job 61635379 has been allocated resources
salloc.exe: Granted job allocation 61635379
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3136 are ready for job
srun: error: x11: no local DISPLAY defined, skipping

[user@cn3136 ~]$ **module load qtltools**
[+] Loading qtltools  1.2  on cn3136
[+] Loading gcc  9.2.0  ...
[+] Loading GSL 2.6 for GCC 9.2.0 ...
[-] Unloading gcc  9.2.0  ...
[+] Loading gcc  9.2.0  ...
[+] Loading openmpi 3.1.4  for GCC 9.2.0
[+] Loading ImageMagick  7.0.8  on cn3136
[+] Loading HDF5  1.10.4
[-] Unloading gcc  9.2.0  ...
[+] Loading gcc  9.2.0  ...
[+] Loading NetCDF 4.7.4_gcc9.2.0
[+] Loading pandoc  2.10  on cn3136
[+] Loading pcre2 10.21  ...
[+] Loading R 4.0.0

[user@cn3136 ~]$ **cd /lscratch/$SLURM\_JOB\_ID**

[user@cn3136 61635379]$ **wget http://jungle.unige.ch/QTLtools\_examples/HG00381.chr22.bam && \
 wget http://jungle.unige.ch/QTLtools\_examples/HG00381.chr22.bam.bai && \
 wget http://jungle.unige.ch/QTLtools\_examples/gencode.v19.annotation.chr22.gtf.gz**
--2020-07-22 15:44:12--  http://jungle.unige.ch/QTLtools_examples/HG00381.chr22.bam
Resolving dtn01-e0 (dtn01-e0)... 10.1.200.237
Connecting to dtn01-e0 (dtn01-e0)|10.1.200.237|:3128... connected.
Proxy request sent, awaiting response... 200 OK
Length: 47167929 (45M) [text/plain]
Saving to: ‘HG00381.chr22.bam’

100%[============================================================>] 47,167,929  16.2MB/s   in 2.8s

2020-07-22 15:44:15 (16.2 MB/s) - ‘HG00381.chr22.bam’ saved [47167929/47167929]

--2020-07-22 15:44:15--  http://jungle.unige.ch/QTLtools_examples/HG00381.chr22.bam.bai
Resolving dtn01-e0 (dtn01-e0)... 10.1.200.237
Connecting to dtn01-e0 (dtn01-e0)|10.1.200.237|:3128... connected.
Proxy request sent, awaiting response... 200 OK
Length: 50136 (49K) [text/plain]
Saving to: ‘HG00381.chr22.bam.bai’

100%[============================================================>] 50,136       173KB/s   in 0.3s

2020-07-22 15:44:15 (173 KB/s) - ‘HG00381.chr22.bam.bai’ saved [50136/50136]

--2020-07-22 15:44:15--  http://jungle.unige.ch/QTLtools_examples/gencode.v19.annotation.chr22.gtf.gz
Resolving dtn01-e0 (dtn01-e0)... 10.1.200.237
Connecting to dtn01-e0 (dtn01-e0)|10.1.200.237|:3128... connected.
Proxy request sent, awaiting response... 200 OK
Length: 881813 (861K) [application/x-gzip]
Saving to: ‘gencode.v19.annotation.chr22.gtf.gz’

100%[============================================================>] 881,813      194KB/s   in 4.4s

2020-07-22 15:44:20 (194 KB/s) - ‘gencode.v19.annotation.chr22.gtf.gz’ saved [881813/881813]


[user@cn3136 61635379]$ **QTLtools quan \
 --bam HG00381.chr22.bam \
 --gtf gencode.v19.annotation.chr22.gtf.gz \
 --sample HG00381 \
 --out-prefix HG00381 \
 --filter-mapping-quality 150 \
 --filter-mismatch 5 \
 --filter-mismatch-total 5 \
 --rpkm**

QTLtools
  * Authors : Olivier DELANEAU / Halit ONGEN / Emmanouil DERMITZAKIS
  * Contact : olivier.delaneau@gmail.com / halit.ongen@unige.ch / Emmanouil.Dermitzakis@unige.ch
  * Webpage : https://qtltools.github.io/qtltools/
  * Version : 1.2
  * Date    : 22/07/2020 - 15:44:33
  * Citation: A complete tool set for molecular QTL discovery and analysis, https://doi.org/10.1038/ncomms15452

QUANTIFY GENES AND EXONS FROM BAM FILES

WARNING: OUTPUT IS NOT COMPATABLE WITH QUANTIFICATIONS GENERATED BEFORE VERSION 1.2

  * Minimum mapping quality: 150
  * Maximum mismatch count per mate-pair: 5
  * Maximum mismatch count per read: 5
  * Not checking properly paired flag
  * Not checking if all blocks of a split read are consistent with the annotation
  * Not filtering reads flagged as duplicate
  * Not filtering reads flagged as failing QC
  * Merging overlapping mate pairs
  * Excluding exons smaller than 0 bp only in exon quantifications
  * Including all gene types
  * Unique hash for this combination of options and GTF file: 2bdBW71i0Bq

Initialize random number generator
  * Built-in seed is 15112011
  * First Integer = 30554
  * First Double = 0.0126695

Reading exons in [gencode.v19.annotation.chr22.gtf.gz]

Opening BAM file [HG00381.chr22.bam]
  * reading header
  * reading index file
  * reading BAM file
  * Expecting 610860 lines
[======================================================================] 100 %
  * DONE: 610860 lines read

WARNING: 8842 unmatched mate pairs found

Printing counts

Printing RPKM

Printing stats

Running time: 3 seconds

[user@cn3136 61635379]$ **exit**
exit
salloc.exe: Relinquishing job allocation 61635379
salloc.exe: Job allocation 61635379 has been revoked.

[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. qtltools.sh). For example:



```

#!/bin/bash

module load qtltools
cd /lscratch/$SLURM_JOB_ID
wget http://jungle.unige.ch/QTLtools_examples/HG00381.chr22.bam && \
    wget http://jungle.unige.ch/QTLtools_examples/HG00381.chr22.bam.bai && \
    wget http://jungle.unige.ch/QTLtools_examples/gencode.v19.annotation.chr22.gtf.gz
QTLtools quan \
    --bam HG00381.chr22.bam \
    --gtf gencode.v19.annotation.chr22.gtf.gz \
    --sample HG00381 \
    --out-prefix HG00381 \
    --filter-mapping-quality 150 \
    --filter-mismatch 5 \
    --filter-mismatch-total 5 \
    --rpkm

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] [--gres=lscratch:#] qtltools.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. qtltools.swarm). For example:



```

QTLtools mbv --bam samp1.bam --vcf samp1.vcf.gz --filter-mapping-quality 150 --out samp1.bamstat.txt
QTLtools mbv --bam samp2.bam --vcf samp2.vcf.gz --filter-mapping-quality 150 --out samp2.bamstat.txt
QTLtools mbv --bam samp3.bam --vcf samp3.vcf.gz --filter-mapping-quality 150 --out samp3.bamstat.txt
QTLtools mbv --bam samp4.bam --vcf samp4.vcf.gz --filter-mapping-quality 150 --out samp4.bamstat.txt

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f qtltools.swarm [-g #] [-t #] --module qtltools
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module qtltools Loads the qtltools module for each subjob in the swarm 
 | |
 | |
 | |








