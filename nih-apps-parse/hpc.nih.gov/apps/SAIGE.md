

document.querySelector('title').textContent = 'SAIGE on Biowulf';
SAIGE on Biowulf


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


SAIGE is installed as a container with it's own R environment on the Biowulf Cluster, please do not load R module when running SAIGE.
If there are conflicts/errors about R, please check the loaded modules with 'module list'.




SAIGE is an R package developed with Rcpp for genome-wide association tests in large-scale data sets and biobanks. The method:



* accounts for sample relatedness based on the generalized mixed models
* allows for model fitting with either full or sparse genetic relationship matrix (GRM)
* works for quantitative and binary traits
* handles case-control imbalance of binary traits
* computationally efficient for large data sets
* performs single-variant association tests
* provides effect size estimation through Firth's Bias-Reduced Logistic Regression
* performs conditional association analysis



SAIGE-GENE (now known as SAIGE-GENE+) are new method extension in the R package for testing rare variant in set-based tests.



* performs BURDEN, SKAT, and SKAT-O tests
* allows for tests on multiple minor allele frequencies cutoffs and functional annotations
* allows for specifying weights for markers in the set-based tests
* performs conditional analysis to identify associations independent from nearly GWAS signals



The package takes genotype file input in the following formats



* PLINK (bed, bim, fam), BGEN, VCF, BCF, SAV


### References:


* Zhou W et.al. *Efficiently controlling for case-control imbalance and sample relatedness in large-scale genetic association studies.*
 Nat Genet. 2018 Sep;50(9):1335-1341.
 [PubMed](https://www.ncbi.nlm.nih.gov/pubmed/30104761) | 
 [Journal](https://www.nature.com/articles/s41588-018-0184-y)


Documentation
* SAIGE Main Site:[Main Site](https://github.com/saigegit/SAIGE)
* SAIGE documentation: [SAIGE-doc](https://saigegit.github.io//SAIGE-doc/)


Important Notes
* Module Name: SAIGE (see [the modules page](/apps/modules.html) for more information)
 
```

step1_fitNULLGLMM.R --help
	

```
* Environment variables set 
	+ $SAIGE\_TEST\_DATA



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive --cpus-per-task=6 --mem=4G**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load SAIGE**
[user@cn3144 ~]$ **cp -r ${SAIGE\_TEST\_DATA:-none}/extdata .**
[user@cn3144 ~]$ **cd extdata**
[user@cn3144 ~]$ **step1\_fitNULLGLMM.R \
 --plinkFile=./input/nfam\_100\_nindep\_0\_step1\_includeMoreRareVariants\_poly \
 --phenoFile=./input/pheno\_1000samples.txt\_withdosages\_withBothTraitTypes.txt \
 --phenoCol=y\_binary \
 --covarColList=a9 \
 --sampleIDColinphenoFile=IID \
 --traitType=binary \
 --outputPrefix=./output/example\_binary\_includenonAutoforvarRatio \
 --nThreads=4 \
 --LOCO=FALSE \
 --relatednessCutoff=0.0 \
 --FemaleCode=2 \
 --MaleCode=1 \
 --IsOverwriteVarianceRatioFile=TRUE** 
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. SAIGE.sh). For example:



```

#!/bin/bash
#SBATCH --cpus-per-task=6
#SBATCH --mem=4G
#SBATCH --time=2:00:00
#SBATCH --partition=norm

set -e
module load SAIGE
step1_fitNULLGLMM.R     \
        --plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly \
        --phenoFile=./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt \
        --phenoCol=y_binary \
        --covarColList=a9 \
        --sampleIDColinphenoFile=IID \
        --traitType=binary        \
        --outputPrefix=./output/example_binary_includenonAutoforvarRatio \
        --nThreads=4    \
        --LOCO=FALSE    \
        --relatednessCutoff=0.0 \
        --FemaleCode=2 \
        --MaleCode=1 \
        --IsOverwriteVarianceRatioFile=TRUE


```

 Submit the job:

```
sbatch SAIGE.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. job.swarm). For example:



```


cd dir1; step1_fitNULLGLMM.R --help 
cd dir2; step1_fitNULLGLMM.R --help

    
```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f job.swarm [-g #] --module SAIGE
```

where
 

|  |  |
| --- | --- |
| -g *#*  | Number of Gigabytes of memory required for each process
 (1 line in the swarm command file)  |
| --module  | Loads the module for each subjob in the swarm  |










