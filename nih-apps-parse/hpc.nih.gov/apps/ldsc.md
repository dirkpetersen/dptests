

document.querySelector('title').textContent = 'LDSC on Biowulf';
LDSC on Biowulf


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



LDSC is a command line tool for estimating heritability and genetic correlation from GWAS summary statistics. LDSC also computes LD Scores.



### References:



The following is from <https://github.com/bulik/ldsc#citation>

* If you use the software or the LD Score regression intercept, please cite

Bulik-Sullivan, et al. LD Score Regression Distinguishes Confounding from Polygenicity in Genome-Wide Association Studies. Nature Genetics, 2015.
* For genetic correlation, please also cite

Bulik-Sullivan, et al. An Atlas of Genetic Correlations across Human Diseases and Traits. bioRxiv doi: http://dx.doi.org/10.1101/014498
* For partitioned heritability, please also cite

Finucane, HK, et al. Partitioning Heritability by Functional Category using GWAS Summary Statistics. bioRxiv doi: http://dx.doi.org/10.1101/014241
* If you find the fact that LD Score regression approximates HE regression to be conceptually useful, please cite

Bulik-Sullivan, Brendan. Relationship between LD Score and Haseman-Elston, bioRxiv doi http://dx.doi.org/10.1101/018283
* For LD Hub, please cite

Zheng, et al. LD Hub: a centralized database and web interface to perform LD score regression that maximizes the potential of summary level GWAS data for SNP heritability and genetic correlation analysis. Bioinformatics (2016) https://doi.org/10.1093/bioinformatics/btw613




Documentation
* [LDSC Home Page](https://github.com/bulik/ldsc)
* [LDSC Documentation](https://github.com/bulik/ldsc/wiki)


Important Notes
* Module Name: ldsc (see [the modules page](/apps/modules.html) for more information)




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

[user@cn3144 ~]$ **module load ldsc**
[+] Loading python 2.7 ...
[+] Loading ldsc, version 1.0.0-92-gcf1707e...
[user@cn3144 ~]$ **wget https://data.broadinstitute.org/alkesgroup/LDSCORE/1kg\_eur.tar.bz2**
--2017-08-01 17:06:24--  https://data.broadinstitute.org/alkesgroup/LDSCORE/1kg_eur.tar.bz2
Resolving dtn01-e0... 10.1.200.237
Connecting to dtn01-e0|10.1.200.237|:3128... connected.
Proxy request sent, awaiting response... 200 OK
Length: 1134910 (1.1M) [application/x-bzip2]
Saving to: “1kg_eur.tar.bz2”

100%[==========================================================================================================>] 1,134,910    833K/s   in 1.3s    

2017-08-01 17:06:30 (833 KB/s) - “1kg_eur.tar.bz2” saved [1134910/1134910]

[user@cn3144 ~]$ **tar -xf 1kg\_eur.tar.bz2**
[user@cn3144 ~]$ **cd 1kg\_eur/**
[user@cn3144 ~]$ **ldsc.py --bfile 22 --l2 --ld-wind-cm 1 --out 22**
*********************************************************************
* LD Score Regression (LDSC)
* Version 1.0.0
* (C) 2014-2015 Brendan Bulik-Sullivan and Hilary Finucane
* Broad Institute of MIT and Harvard / MIT Department of Mathematics
* GNU General Public License v3
*********************************************************************
Call: 
./ldsc.py \
--ld-wind-cm 1.0 \
--out 22 \
--bfile 22 \
--l2  

Beginning analysis at Tue Aug  1 17:12:04 2017
Read list of 19156 SNPs from 22.bim
Read list of 379 individuals from 22.fam
Reading genotypes from 22.bed
After filtering, 19156 SNPs remain
Estimating LD Score.
Writing LD Scores for 19156 SNPs to 22.l2.ldscore.gz

Summary of LD Scores in 22.l2.ldscore.gz
         MAF        L2
mean  0.2323   18.5353
std   0.1453   16.1039
min   0.0013    0.0657
25%   0.1042    7.8392
50%   0.2243   13.4837
75%   0.3549   22.9722
max   0.5000  109.7163

MAF/LD Score Correlation Matrix
        MAF      L2
MAF  1.0000  0.2749
L2   0.2749  1.0000
Analysis finished at Tue Aug  1 17:12:08 2017
Total time elapsed: 3.31s
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. ldsc.sh). For example:



```

#!/bin/sh
module load ldsc

wget https://data.broadinstitute.org/alkesgroup/LDSCORE/1kg_eur.tar.bz2
tar -xf 1kg_eur.tar.bz2
cd 1kg_eur/
ldsc.py --bfile 22 --l2 --ld-wind-cm 1 --out 22

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] ldsc.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit set of independent commands requiring identical resources.
Create a swarmfile (e.g. ldsc.swarm). For example:



```

cd /data/$USER/ldsc/set1 && ldsc.py --bfile 22 --l2 --ld-wind-cm 1 --out 22
cd /data/$USER/ldsc/set2 && ldsc.py --bfile 22 --l2 --ld-wind-cm 1 --out 22
cd /data/$USER/ldsc/set3 && ldsc.py --bfile 22 --l2 --ld-wind-cm 1 --out 22
cd /data/$USER/ldsc/set4 && ldsc.py --bfile 22 --l2 --ld-wind-cm 1 --out 22

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f ldsc.swarm [-g #] [-t #] --module ldsc
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g #  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t #  Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module ldsc  Loads the TEMPLATE module for each subjob in the swarm 
 | |
 | |
 | |










