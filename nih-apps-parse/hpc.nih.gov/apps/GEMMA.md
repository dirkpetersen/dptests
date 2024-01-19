

document.querySelector('title').textContent = ' GEMMA on Biowulf & Helix';
 GEMMA on Biowulf & Helix


|  |
| --- |
| 
Quick Links
[Documentation & Notes](#doc)
[Interactive job](#int-threaded)
[Batch job](#sbatch-threaded)
 |



GEMMA is the software implementing the Genome-wide Efficient Mixed Model Association algorithm for a standard linear mixed model and some of its close relatives for genome-wide association studies (GWAS):
* It fits a univariate linear mixed model (LMM) for marker association tests with a single phenotype to account for population stratification and sample structure, and for estimating the proportion of variance in phenotypes explained (PVE) by typed genotypes (i.e. "chip heritability").
* It fits a multivariate linear mixed model (mvLMM) for testing marker associations with multiple phenotypes simultaneously while controlling for population stratification, and for estimating genetic correlations among complex phenotypes.
* It fits a Bayesian sparse linear mixed model (BSLMM) using Markov chain Monte Carlo (MCMC) for estimating PVE by typed genotypes, predicting phenotypes, and identifying associated markers by jointly modeling all markers while controlling for population structure.
* It estimates variance component/chip heritability, and partitions it by different SNP functional categories. In particular, it uses HE regression or REML AI algorithm to estimate variance components when individual-level data are available. It uses MQS to estimate variance components when only summary statisics are available.


GEMMA was developed in the Zhou Lab at U. Michigan. GEMMA website.



Documentation and Notes
* For Documentation, see [GEMMA Github](https://github.com/genetics-statistics/GEMMA)
* Run `$ gemma -h` to see commandline help documentations
* Module Name: GEMMA (see [the modules page](/apps/modules.html) for more information)
* Use environment variable $GEMMA\_TEMPLATES to access GEMMA example files



Interactive job 
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf ~]$ **sinteractive**
salloc.exe: Pending job allocation 57710203
salloc.exe: job 57710203 queued and waiting for resources
salloc.exe: job 57710203 has been allocated resources
salloc.exe: Granted job allocation 57710203
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3115 are ready for job

[user@cn3115 ~]$ **mkdir -p /data/$USER/gemma-test**

[user@cn3115 ~]$ **cd /data/$USER/gemma-test**

[user@cn3115 ~]$ **module load GEMMA**
[+] Loading GEMMA 0.98.1  ...

[user@cn3115 ~]$ **gemma -g $GEMMA\_EXAMPLES/mouse\_hs1940.geno.txt.gz -p $GEMMA\_EXAMPLES/mouse\_hs1940.pheno.txt -gk -o mouse\_hs1940**
GEMMA 0.98.1 (2018-12-10) by Xiang Zhou and team (C) 2012-2018
Reading Files ...
## number of total individuals = 1940
## number of analyzed individuals = 1410
## number of covariates = 1
## number of phenotypes = 1
## number of total SNPs/var        =    12226
## number of analyzed SNPs         =    10768
Calculating Relatedness Matrix ...
================================================== 100%
**** INFO: Done.

[user@cn3115 ~]$ **exit**
exit
salloc.exe: Relinquishing job allocation 57710203
salloc.exe: Job allocation 57710203 has been revoked.

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. gemma.sh). For example, to run the threaded version of hyphy:



```

#!/bin/bash
# this file is called gemma.sh

module load GEMMA
mkdir -p /data/$USER/gemma-test
cd /data/$USER/gemma-test

gemma -g $GEMMA_EXAMPLES/mouse_hs1940.geno.txt.gz \
      -p $GEMMA_EXAMPLES/mouse_hs1940.pheno.txt \
      -gk -o mouse_hs1940

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
$ sbatch --cpus-per-task=# [--mem=#] gemma.sh
```







