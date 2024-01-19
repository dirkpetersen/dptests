

document.querySelector('title').textContent = 'Eagle: reference-based haplotype phasing';
**Eagle: reference-based haplotype phasing**


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



Eagle performs a reference-based haplotype phasing. It attains high
accuracy across a broad range of cohort sizes by efficiently leveraging information 
from large external reference panels (such as the Haplotype Reference Consortium; HRC) 
using a new data structure based on the positional Burrows-Wheeler transform.



### References:


* Po-Ru Loh, Petr Danecek, Pier Francesco Palamara, Christian Fuchsberger, Yakir A Reshef,
Hilary K Finucane, Sebastian Schoenherr, Lukas Forer, Shane McCarthy, Goncalo R Abecasis,
Richard Durbin and Alkes L Price, 
"Reference-based phasing using the Haplotype Reference Consortium panel", Nature Genetics, 2016, **48**(11), 1443-1450.

* Po-Ru Loh, Pier Francesco Palamara and Alkes L Price, "Fast and accurate long-range phasing in a UK Biobank cohort", Nature Genetics, 2016, **48**(7), 811-819.


Documentation
* [Eagle Home Page](https://data.broadinstitute.org/alkesgroup/Eagle)
* [Eagle GitHub Page](http://github.com/poruloh/Eagle)


Important Notes
* Module Name: Eagle (see [the modules page](/apps/modules.html) for more information)
* Unusual environment variables set
	+ **EAGLE\_DIR** Eagle installation directory
	+ **EAGLE\_BIN** Eagle executable folder
	+ **EAGLE\_DATA** sample data for running Eagle
	+ **EAGLE\_TABLES** sample data for running Eagle



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive**
[user@cn3144 ~]$ **module load Eagle** 

```

Copy sample data from an application folder to current folder:

```

[user@cn3144 ~]$ **cp $EAGLE\_DATA/\* .** 

```

Run Eagle on the sample data: :

```

[user@cn3144 ~]$ **eagle --bfile=EUR\_test --geneticMapFile=USE\_BIM --chrom=21 --outPrefix=phased --numThreads=4 2>&1 | tee example.log** 

                      +-----------------------------+
                      |                             |
                      |   Eagle v2.4                |
                      |   December 13, 2017         |
                      |   Po-Ru Loh                 |
                      |                             |
                      +-----------------------------+

Copyright (C) 2015-2017 Harvard University.
Distributed under the GNU GPLv3+ open source license.

Command line options:

eagle \
    --bfile=EUR_test \
    --geneticMapFile=USE_BIM \
    --chrom=21 \
    --outPrefix=phased \
    --numThreads=4

Setting number of threads to 4

=== Reading genotype data ===

Reading fam file: EUR_test.fam
Total indivs in PLINK data: Nbed = 379
Total indivs stored in memory: NpreQC = 379
Reading bim file: EUR_test.bim
Total snps in PLINK data: Mbed = 2000
Restricting to 1813 SNPs on chrom 21 in region [bpStart,bpEnd] = [0,1e+09]
Total SNPs stored in memory: MpreQC = 1813
Allocating 1813 x 379 bytes to temporarily store genotypes
Reading genotypes and performing QC filtering on snps and indivs...
Reading bed file: EUR_test.bed
    Expecting 190000 (+3) bytes for 379 indivs, 2000 snps

Total post-QC indivs: N = 379
Total post-QC SNPs: M = 1813
MAF spectrum:
     0- 5%:     495
     5-10%:     290
    10-20%:     332
    20-30%:     248
    30-40%:     234
    40-50%:     214
Physical distance range: 9752235 base pairs
Genetic distance range:  23.0881 cM
Average # SNPs per cM:   79
Auto-selecting --maxBlockLen: 0.25 cM
Number of <=(64-SNP, 0.25cM) segments: 68
Average # SNPs per segment: 26
Estimating LD scores using 379 indivs
Fraction of heterozygous genotypes: 0.246308
Typical span of default 100-het history length: 5.17 cM
Setting --histFactor=1.00

BEGINNING STEP 1

Time for step 1: 0.867686
Time for step 1 MN^2: 0.0521836

Making hard calls (time: 0.0207999)


BEGINNING STEP 2

BATCH 1 OF 1
Building hash tables
.................................................................. (time: 0.136335)

Phasing samples 1-379
Time for phasing batch: 1.03954

Making hard calls (time: 0.020123)

Time for step 2: 1.19602
Time for step 2 MN^2: 0.158607


BEGINNING STEP 3 (PBWT ITERS)

Auto-selecting number of PBWT iterations: setting --pbwtIters to 2


BEGINNING PBWT ITER 1

BATCH 1 OF 10

Phasing samples 1-37
Time for phasing batch: 3.31806

BATCH 2 OF 10

Phasing samples 38-75
Time for phasing batch: 3.23385

...

BATCH 10 OF 10

Phasing samples 342-379
Time for phasing batch: 3.21097

Time for PBWT iter 1: 31.8771

BEGINNING PBWT ITER 2

BATCH 1 OF 10

Phasing samples 1-37
Time for phasing batch: 5.23776

BATCH 2 OF 10

Phasing samples 38-75
Time for phasing batch: 5.15485

...

BATCH 9 OF 10

Phasing samples 304-341
Time for phasing batch: 5.06871

BATCH 10 OF 10

Phasing samples 342-379
Time for phasing batch: 5.19495

Time for PBWT iter 2: 51.1316
Writing .haps.gz and .sample output
Time for writing output: 0.23035
Total elapsed time for analysis = 85.4332 sec

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. Eagle.sh). For example:



```

#!/bin/bash
module load Eagle     
eagle \
    --vcf=EUR_test.vcf.gz \
    --geneticMapFile=$EAGLE_TABLES/genetic_map_hg19_withX.txt.gz \
    --chrom=21 \
    --outPrefix=phased \
    --numThreads=4 \
    2>&1 | tee example_vcf.log

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] Eagle.sh
```







