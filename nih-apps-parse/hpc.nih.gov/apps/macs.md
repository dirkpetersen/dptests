

document.querySelector('title').textContent = 'macs on Biowulf';
macs on Biowulf


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


Model-based Analysis of ChIP-Seq (MACS) is used on short reads sequencers
such as Genome Analyzer (Illumina / Solexa). MACS empirically models the length
of the sequenced ChIP fragments, which tends to be shorter than sonication or
library construction size estimates, and uses it to improve the spatial
resolution of predicted binding sites. MACS also uses a dynamic Poisson
distribution to effectively capture local biases in the genome sequence,
allowing for more sensitive and robust prediction. MACS compares favorably to
existing ChIP-Seq peak-finding algorithms and can be used for ChIP-Seq with or
without control samples.


### References:


* Y. Zhang, T. Liu, C. A. Meyer, J. Eeckhoute, D. S. Johnson, B. E. Bernstein, C. Nusbaum, R. M. Myers, M. Brown, W. Li and X. S. Liu. *Model-based Analysis of ChIP-Seq (MACS)* Genome Biology 2008. 
[PubMed](https://www.ncbi.nlm.nih.gov/pubmed/18798982) | 
[PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2592715/) | 
[Journal](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2008-9-9-r137)


Documentation
* Source code repository: [on GitHub](https://github.com/taoliu/MACS)
* Manual [on GitHub](https://github.com/taoliu/MACS/wiki)


Important Notes
* Module Name: macs (see [the modules page](/apps/modules.html) for more information)
* Example data in `$MACS_TEST_DATA`
* MACS will need lscratch for storage temporary files.



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=4g --gres=lscratch:5** 
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load macs**
[+] Loading macs  2.2.6

[user@cn3144 ~]$ **macs2**
usage: macs2 [-h] [--version]
             {callpeak,bdgpeakcall,bdgbroadcall,bdgcmp,bdgopt,cmbreps,bdgdiff,filterdup,predictd,pileup,randsample,refinepeak}
             ...

[user@cn3144 ~]$ **cp $MACS\_TEST\_DATA/\*.bam .**
[user@cn3144 ~]$ **ls -lh**
-rw-rw-r-- 1 user group 387M Feb 14 10:19 ENCFF001NGB.bam
-rw-rw-r-- 1 user group 458M Feb 14 10:19 ENCFF001NHS.bam

[user@cn3144 ~]$ # ENCFF001NHS.bam is the control (IgG) data
[user@cn3144 ~]$ **macs2 callpeak \
 -t ENCFF001NGB.bam -c ENCFF001NHS.bam \
 -n test -g mm**

INFO  @ Wed, 14 Feb 2018 09:47:20:
# Command line: callpeak -t /usr/local/apps/macs/TEST_DATA/ENCFF001NGB.bam -c /usr/local/apps/macs/TEST_DATA/ENCFF001NH
S.bam -n test -g mm
# ARGUMENTS LIST:
# name = test
# format = AUTO
# ChIP-seq file = ['/usr/local/apps/macs/TEST_DATA/ENCFF001NGB.bam']
# control file = ['/usr/local/apps/macs/TEST_DATA/ENCFF001NHS.bam']
# effective genome size = 1.87e+09
# band width = 300
# model fold = [5, 50]
# qvalue cutoff = 5.00e-02
# Larger dataset will be scaled towards smaller dataset.
# Range for calculating regional lambda is: 1000 bps and 10000 bps
# Broad region calling is off
# Paired-End mode is off

INFO  @ Wed, 14 Feb 2018 09:47:20: #1 read tag files...
[...snip...]

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. macs.sh), which uses the input file 'macs.in'. For example:



```

#! /bin/bash

module load macs || exit 1
macs2 callpeak \
    -t treatment1.bam \
    -c control.bam \
    --call-summits \
    --tempdir /lscratch/$SLURM_JOB_ID \
    --qvalue 0.01 \
    -n treatment1 -g mm

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=2 --mem=10g  --gres=lscratch:20 macs.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. macs.swarm). For example:



```

macs2 callpeak -t treatment1.bam -c control.bam --name treatment1
macs2 callpeak -t treatment2.bam -c control.bam --name treatment2
macs2 callpeak -t treatment3.bam -c control.bam --name treatment3

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f macs.swarm -g 10 --module macs --gres=lscratch:10
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g #  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t #  Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module macs  Loads the macs module for each subjob in the swarm 
 | |
 | |
 | |








