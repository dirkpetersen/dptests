

document.querySelector('title').textContent = 'mach2qtl on Biowulf';
mach2qtl on Biowulf



|  |
| --- |
| 
Quick Links
[Interactive mach2qtl on Biowulf](#int)
[Serial mach2qtl on Biowulf](#serial)

[Swarm of mach2qtl jobs](#swarm)
[Documentation](#doc)
 |





Description
Mach2qtl uses dosages/posterior probabilities inferred with
[mach](https://hpc.cit.nih.gov/apps/Mach1.0.html) or
[minimac](https://hpc.cit.nih.gov/apps/minimac.html) as predictors in a linear
regression to test association with quantitative traits.


### References


* Yun Li, C. J. Willer, J. Ding, P. Scheet and G. R. Abecasis. *MaCH: using 
 sequence and genotype data to estimate haplotypes and unobserved genotypes*. 
 Genet Epidemiol 2010, 34:816-834.
 [Pubmed](http://www.ncbi.nlm.nih.gov/pubmed/21058334) | 
 [PMC](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3175618/) | 
 [Journal](http://onlinelibrary.wiley.com/doi/10.1002/gepi.20533/abstract)
* Yu-Fang Pei, L. Zhang, J. Li, H. Deng. *Analyses and Comparison of 
 Imputation-Based Association Methods*. PLoS One 2010, 5:e10827.
 [Pubmed](http://www.ncbi.nlm.nih.gov/pubmed/20520814) | 
 [PMC](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2877082/) | 
 [Journal](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0010827)


### Web sites


* [Home page](http://www.unc.edu/~yunmli/software.html)
* [Previous home page](http://csg.sph.umich.edu/abecasis/MACH/download/)



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

[user@cn3144 ~]$ **module load mach2qtl**

[user@cn3144 ~]$ **cd /data/$USER/test\_data/mach2qtl**

[user@cn3144 ~]$ **mach2qtl -d sample.dat -p sample.ped -i sample.mlinfo \
 --probfile sample.mlprob \
 --samplesize > test.out**

[user@cn3144 ~]$ **less test.out**
Mach2Qtl V1.1.3 (2013-04-25) -- QTL Association Mapping with Imputed Allele Counts
(c) 2007 Goncalo Abecasis, Yun Li


The following parameters are in effect:

Available Options
         Phenotypic Data : --datfile [sample.dat], --pedfile [sample.ped]
   Imputed Allele Counts : --infofile [sample.mlinfo], --dosefile [],
                           --probfile [sample.mlprob]
        Analysis Options : --useCovariates [ON], --quantileNormalization,
                           --dominant, --recessive, --additive [ON]
                  Output : --samplesize [ON], --rsqcutoff [0.00]

FITTED MODELS (for covariate adjusted residuals)
======================================================
          Trait     Raw Mean Raw Variance         Mean     Variance
        lg10CRP      0.52383      0.13666     -0.00000      0.13180

Loading marker information ...
   1001 markers will be analyzed

Processing prob file ...

Executing first [ass analysis of imputed genotypes ...
Executing second pass analysis of imputed genotypes ...

                    INFORMATION FROM .info FILE             QTL ASSOCIATION ADDITIVE
                    =====================================   =================================
TRAIT               MARKER               ALLELES  FREQ1  RSQR   EFFECT2  STDERR  CHISQ     PVALUE N
lg10CRP             rs9977217            A,G      .3950  .9805   0.020   0.018   1.2320     0.267    867
lg10CRP             rs9977094            A,C      .9894  .3454   0.179   0.133   1.8130    0.1782    867
lg10CRP             rs10854272           C,T      .4169  .9562   0.014   0.018   0.6349    0.4256    867
...
Analysis took 1 seconds

[user@cn3144 ~]$ **mach2qtl -d sample.dat -p sample.ped -i sample.mlinfo \
 --probfile sample.mlprob --dominant --recessive \
 --samplesize > test.out**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Running a single mach2qtl batch job on Biowulf


To run a single `mach2qtl` analysis on biowulf, a batch script similar to the
exammple below will be required:



```

#! /bin/bash
set -e

cd /data/$USER/test_data/mach2qtl
module load mach2qtl/1.1.3
mach2qtl -d sample.dat -p sample.ped -i sample.mlinfo \
    --probfile sample.mlprob --dominant --recessive \
    --samplesize > test.out

```

This script can then be submitted to the slurm batch queue with a command like



```

biowulf$ **sbatch --mem=2g batch\_script.sh**

```


Running a swarm of mach2qtl batch jobs on Biowulf


It is easiest to run a large number of simultaneous `mach2qtl` jobs via swarm.
Set up a swarm command file along the following lines:



```

#------- this file is swarmcmd ------------------
mach2qtl -d sample1.dat -p sample1.ped -i sample1.mlinf --probfile sample1.mlp  > sample1.out
mach2qtl -d sample2.dat -p sample2.ped -i sample2.mlinf --probfile sample2.mlp  > sample2.out
mach2qtl -d sample1.dat -p sample3.ped -i sample3.mlinf --probfile sample3.mlp  > sample3.out
mach2qtl -d sample4.dat -p sample4.ped -i sample4.mlinf --probfile sample4.mlp  > sample4.out
[...]

```

If each `mach2qtl` process requires less than 1 GB of memory, submit this to the batch system with the command:



```

swarm -f swarmcmd --module mach2qtl

```

If each `mach2qtl` process requires more than 1 GB of memory, use



```

swarm -g # -f cmdfile --module mach2qtl

```

where '#' is the number of Gigabytes of memory required by each `mach2qtl`
process. The swarm program will package the commands for best efficiency and
send them to the batch system as a job array.



Documentation


Mach2qtl readme:



```


mach2qtl
========
QTL analysis based on imputed dosages/posterior_probabilities


Options:
--------

    --datfile : merlin marker info file for phenotypes (required)

    --pedfile : merlin pedigree file for phenotypes (required)

    --infofile : mach1 output .info or .mlinfo file (required)

    --dosefile : mach1 output .dose or .mldose file (required if probfile not available)

    --probfile : mach1 output .prob or .mlprob file (required for non-additive model)

    --useCovariates : adjust for covariates (default is ON)

    --quantileNormalization : apply inverse normal transformation to quantitative trait(s) (default if OFF)

    --dominant : dominant model (default if OFF)

    --recessive : recessive model (default if OFF)

    --additive : additive model (default if ON)

Sample command line & Example
-----------------------------

    see sub-folder examples/

Coding
======
For recessive model: AL1/AL1                    1
                     AL1/AL2 & AL2/AL2          0
    so a positive beta means AL1/AL1 is associated with higher trait values)

For dominant model:  AL2/AL2                    1
                     AL1/AL2 & AL1/AL1          0
    so a positive beta means AL2/AL2 is associated with higher trait values)

For additive model:  AL1            0
             AL2            1
             (notice the negative sign in "double effect = -numerator / denominator;"
    so a positive beta means AL2 is associated with higher trait values)

Procedure:
----------
Calculate residuals (if related, use VC)
Perform residual-SNP association

Selected Output
----------------
Raw Mean/Variance vs Mean/Variance
        former for y (could be transformed)
        latter for residuals (even if no residuals, Raw Mean and Mean could be different b/c latter would be centered and thus would be ZERO)


ref:
----
Li Y, Willer CJ, Sanna S, and Abecasis GR (2009). Genotype imputation. Annu Rev Genomics Hum Genet. 10: 387-406. 
Chen WM and Abecasis GR (2007). Family-based association tests for genomewide association scans. Am J Hum Genet 81:913-26


```



