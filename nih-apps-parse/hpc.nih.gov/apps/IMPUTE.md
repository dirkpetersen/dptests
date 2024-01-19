

document.querySelector('title').textContent = 'Impute, qctool, gtool and snptest on Biowulf';
Impute, qctool, gtool and snptest on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
[Swarm of jobs](#swarm) 
[Note about parallelization](#par)
 |



**[IMPUTE](https://mathgen.stats.ox.ac.uk/impute/impute.html)** version 2 (also known as IMPUTE2) is a genotype imputation and haplotype phasing program based on ideas from Howie et al. 2009: 

**[QCTOOL](http://www.well.ox.ac.uk/~gav/qctool/#overview)** is a command-line utility program for basic quality control of GWAS datasets. 
It supports the same file formats used by the WTCCC studies, as well as the binary file format described on the qctool webpage and the Variant Call Format, and is 
designed to work seamlessly with SNPTEST and related tools 

**[GTOOL](http://www.well.ox.ac.uk/~cfreeman/software/gwas/gtool.html)** is a program for transforming sets of genotype data for use with 
the programs SNPTEST and IMPUTE. GTOOL can be used to:
* generate subsets of genotype data,
 * to convert genotype data between the PED file format and the file format used by SNPTEST and IMPUTE,
 * merge genotype datasets together,
 * orient genotype data according to a strand file.



**[SNPTEST](https://mathgen.stats.ox.ac.uk/genetics_software/snptest/snptest.html)** is a program for the analysis of single SNP association 
in genome-wide studies. The tests implemented include
* Binary (case-control) phenotypes, single and multiple quantitative phenotypes
 * Bayesian and Frequentist tests
 * Ability to condition upon an arbitrary set of covariates and/or SNPs.
 * Various different methods for the dealing with imputed SNPs.


Impute, qctool, gtool and snptest were all developed at the University of Oxford. 




Documentation
[Impute2 website](https://mathgen.stats.ox.ac.uk/impute/impute_v2.html#Instructions)  

[gtool website](http://www.well.ox.ac.uk/~cfreeman/software/gwas/gtool.html)  

[Snptest website](https://mathgen.stats.ox.ac.uk/genetics_software/snptest/snptest.html)  

[qctool website](http://www.well.ox.ac.uk/~gav/qctool/#overview)

Important Notes
* Module Names: impute, qctool, gtool, snptest (see [the modules page](/apps/modules.html) for more information)
* To add all these modules to your path, use 

```

module load impute qctool gtool snptest  

```



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

[user@cn3144 ~]$ **module load impute**

[user@cn3144 ~]$ **impute2 -ref\_samp\_out -m ./chr16.map -h ./chr16.haps \
 -l ./chr16.legend -g ./chr16.reference.gtypes -s ./chr16.reference.strand \
 -Ne 11418 -int 5000000 5500000 -buffer 250 -k 10 -iter 10 -burnin 3 \
 -o ./Results/chr16.multi\_panel.ref\_gtypes.impute2 \
 -i ./Results/chr16.multi\_panel.ref\_gtypes.impute2.info \
 -r ./Results/chr16.multi\_panel.ref\_gtypes.impute2.summary**

The seed for the random number generator is 1115038504.

Command-line input: impute2 -ref_samp_out -m ./chr16.map -h ./chr16.haps -l ./chr16.legend -g ./chr16.reference.gtypes -s ./chr16.reference.strand -Ne 11418 -int 5000000 5500000 -buffer 250 -k 10 -iter 10 -burnin 3 -o ./Results/chr16.multi_panel.ref_gtypes.impute2 -i ./Results/chr16.multi_panel.ref_gtypes.impute2.info -r ./Results/chr16.multi_panel.ref_gtypes.impute2.summary


======================
 IMPUTE version 2.0.3 
======================

Copyright 2008 Bryan Howie, Peter Donnelly, and Jonathan Marchini
Please see the LICENCE file included with this program for conditions of use.

    haplotypes file : ./chr16.haps
        legend file : ./chr16.legend
 ref genotypes file : NULL
ref gen strand file : NULL
     genotypes file : ./chr16.reference.gtypes
        strand file : ./chr16.reference.strand
           map file : ./chr16.map
 excluded SNPs file : NULL
 included SNPs file : NULL
    ref samp infile : NULL
        output file : ./Results/chr16.multi_panel.ref_gtypes.impute2
          info file : ./Results/chr16.multi_panel.ref_gtypes.impute2.info
       summary file : ./Results/chr16.multi_panel.ref_gtypes.impute2.summary
       [...]
  Accuracy assessment for imputation of type 0 SNPs (those with data in the haploid reference panel only) .The maximum imputed genotype calls are distributed as follows:
  Interval  #Genotypes %Concordance         Interval  %Called %Concordance
  [0.0-0.1]          0          0.0         [ >= 0.0]   100.0         96.7
  [0.1-0.2]          0          0.0         [ >= 0.1]   100.0         96.7
  [0.2-0.3]          0          0.0         [ >= 0.2]   100.0         96.7
  [0.3-0.4]          0          0.0         [ >= 0.3]   100.0         96.7
  [0.4-0.5]         10         40.0         [ >= 0.4]   100.0         96.7

[user@cn3144 ~]$  **ls Results**
chr16.multi_panel.ref_gtypes.impute2
chr16.multi_panel.ref_gtypes.impute2.info
chr16.multi_panel.ref_gtypes.impute2.summary
chr16.multi_panel.ref_gtypes.impute2_refsamp1.gz
chr16.multi_panel.ref_gtypes.impute2_refsamp10.gz
chr16.multi_panel.ref_gtypes.impute2_refsamp2.gz
chr16.multi_panel.ref_gtypes.impute2_refsamp3.gz
chr16.multi_panel.ref_gtypes.impute2_refsamp4.gz
chr16.multi_panel.ref_gtypes.impute2_refsamp5.gz
chr16.multi_panel.ref_gtypes.impute2_refsamp6.gz
chr16.multi_panel.ref_gtypes.impute2_refsamp7.gz
chr16.multi_panel.ref_gtypes.impute2_refsamp8.gz
chr16.multi_panel.ref_gtypes.impute2_refsamp9.gz

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. impute.sh). For example:



```

#!/bin/bash
set -e
module load impute

impute2 -ref_samp_out -m chr16.map -h chr16.haps  -l chr16.legend \
	-g gtypes -s refstrand1  -Ne 11418 -int 5000000 5500000 -buffer 250 \
	-k 10 -iter 10 -burnin 3  -o out1  -i info1  -r summary1

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch  [--mem=#] impute.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. impute.swarm). For example:



```

# this file is impute.swarm
cd /data/user/dir1; impute2 -ref_samp_out -m chr16.map -h chr16.haps  \
      -l chr16.legend -g gtypes -s refstrand1  -Ne 11418 -int 5000000 5500000 \
      -buffer 250 -k 10 -iter 10 -burnin 3  -o out1  -i info1  -r summary1
cd /data/user/dir2; impute2 -ref_samp_out -m chr26.map -h chr26.haps \
      -l chr26.legend -g gtypes -s refstrand2  -Ne 22428 -int 5000000 5500000 \
      -buffer 250 -k 20 -iter 20 -burnin 3  -o out2  -i info2  -r summary2
cd /data/user/dir3; impute2 -ref_samp_out -m chr36.map -h chr36.haps  \
      -l chr36.legend -g gtypes -s refstrand3  -Ne 33438 -int 5000000 5500000 \
      -buffer 250 -k 30 -iter 30 -burnin 3  -o out3  -i info3  -r summary3
[...]

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f impute.swarm [-g #]  --module impute,qctool,gtool
```

where


|  |  |  |  |
| --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | --module impute Loads the impute module for each subjob in the swarm 
 | |
 | |



Note about parallelization

In principle, it is possible to impute genotypes across an entire chromosome in a single run of IMPUTE2. However, we prefer to 
split each chromosome into smaller chunks for analysis, both because the program produces higher accuracy over short genomic 
regions and because imputing a chromosome in chunks is a good computational strategy: the chunks can be imputed in parallel 
on multiple computer processors, thereby decreasing the real computing time and limiting the amount of memory needed for each run.

We therefore recommend using the program on regions of ~5 Mb or shorter, and versions from v2.1.2 onward will throw an error 
if the analysis interval plus buffer region is longer than 7 Mb. People who have good reasons to impute a longer region in a single 
run can override this behavior with the -allow\_large\_regions flag.

See [this informative snippet](https://mathgen.stats.ox.ac.uk/impute/impute_v2.html#whole_chroms) from the Impute website for more details 
about dealing with whole chromosomes.
























