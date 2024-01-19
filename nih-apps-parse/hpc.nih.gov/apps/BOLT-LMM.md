

document.querySelector('title').textContent = 'BOLT-LMM on Biowulf';
BOLT-LMM on Biowulf


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



The BOLT-LMM software package currently consists of two main algorithms, the BOLT-LMM algorithm for mixed model association testing, and the BOLT-REML algorithm for variance components analysis (i.e., partitioning of SNP-heritability and estimation of genetic correlations).



### References:


* Loh, P.-R. et al. Efficient Bayesian mixed model analysis increases association power in large cohorts. Nature Genetics 47, 284–290 (2015).


Documentation
* BOLT-LMM Main Site: <https://data.broadinstitute.org/alkesgroup/BOLT-LMM/>


Important Notes
* Module Name: BOLT-LMM (see [the modules page](/apps/modules.html) for more information)
* Starting with v1.3, multi-threaded support available for BGEN-formatted data files
* Example files in /usr/local/apps/BOLT-LMM/v2.3/example



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session (user input in **bold**):



```

[user@biowulf ~]$ **sinteractive -c4 --mem=8g --gres=lscratch:10**
salloc.exe: Pending job allocation 11198398
salloc.exe: job 11198398 queued and waiting for resources
salloc.exe: job 11198398 has been allocated resources
salloc.exe: Granted job allocation 11198398
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn0852 are ready for job
srun: error: x11: no local DISPLAY defined, skipping
error: unable to open file /tmp/slurm-spank-x11.11198398.0
slurmstepd: error: x11: unable to read DISPLAY value

[user@cn0852 ~]$ **cd /lscratch/$SLURM\_JOB\_ID**

[user@cn0852 11198398]$ **module load BOLT-LMM**
[+] Loading BOLT-LMM  v2.3.4  on cn0852

[user@cn0852 11198398]$ **cp -r /usr/local/apps/BOLT-LMM/v2.3/example .**

[user@cn0852 11198398]$ **cp -r /usr/local/apps/BOLT-LMM/v2.3/tables .**

[user@cn0852 11198398]$ **cd example/**

[user@cn0852 example]$ **bolt \
 --bfile=EUR\_subset \
 --remove=EUR\_subset.remove \
 --exclude=EUR\_subset.exclude \
 --phenoFile=EUR\_subset.pheno.covars \
 --phenoCol=PHENO \
 --covarFile=EUR\_subset.pheno.covars \
 --covarCol=CAT\_COV \
 --qCovarCol=QCOV{1:2} \
 --modelSnps=EUR\_subset.modelSnps \
 --lmm \
 --LDscoresFile=../tables/LDSCORE.1000G\_EUR.tab.gz \
 --numThreads=${SLURM\_CPUS\_PER\_TASK} \
 --statsFile=example.stats \
 --dosageFile=EUR\_subset.dosage.chr17first100 \
 --dosageFile=EUR\_subset.dosage.chr22last100.gz \
 --dosageFidIidFile=EUR\_subset.dosage.indivs \
 --statsFileDosageSnps=example.dosageSnps.stats \
 --impute2FileList=EUR\_subset.impute2FileList.txt \
 --impute2FidIidFile=EUR\_subset.impute2.indivs \
 --statsFileImpute2Snps=example.impute2Snps.stats \
 --dosage2FileList=EUR\_subset.dosage2FileList.txt \
 --statsFileDosage2Snps=example.dosage2Snps.stats \
 2>&1 | tee example.log # log output written to stdout and stderr**
                      +-----------------------------+
                      |                       ___   |
                      |   BOLT-LMM, v2.3.4   /_ /   |
                      |   August 10, 2019     /_/   |
                      |   Po-Ru Loh            //   |
                      |                        /    |
                      +-----------------------------+

Copyright (C) 2014-2019 Harvard University.
Distributed under the GNU GPLv3 open source license.

Compiled with USE_SSE: fast aligned memory access
Compiled with USE_MKL: Intel Math Kernel Library linear algebra
Boost version: 1_58

Command line options:

bolt \
    --bfile=EUR_subset \
    --remove=EUR_subset.remove \
    --exclude=EUR_subset.exclude \
    --phenoFile=EUR_subset.pheno.covars \
    --phenoCol=PHENO \
    --covarFile=EUR_subset.pheno.covars \
    --covarCol=CAT_COV \
    --qCovarCol=QCOV{1:2} \
    --modelSnps=EUR_subset.modelSnps \
    --lmm \
    --LDscoresFile=../tables/LDSCORE.1000G_EUR.tab.gz \
    --numThreads=2 \
    --statsFile=example.stats \
    --dosageFile=EUR_subset.dosage.chr17first100 \
    --dosageFile=EUR_subset.dosage.chr22last100.gz \
    --dosageFidIidFile=EUR_subset.dosage.indivs \
    --statsFileDosageSnps=example.dosageSnps.stats \
    --impute2FileList=EUR_subset.impute2FileList.txt \
    --impute2FidIidFile=EUR_subset.impute2.indivs \
    --statsFileImpute2Snps=example.impute2Snps.stats \
    --dosage2FileList=EUR_subset.dosage2FileList.txt \
    --statsFileDosage2Snps=example.dosage2Snps.stats

Verifying contents of --dosage2FileList: EUR_subset.dosage2FileList.txt
Checking map file EUR_subset.dosage2.chr17first100.map and 2-dosage genotype file EUR_subset.dosage2.chr17first100.gz
Checking map file EUR_subset.dosage2.chr17second100.map and 2-dosage genotype file EUR_subset.dosage2.chr17second100
Checking map file EUR_subset.dosage2.chr22last100.map and 2-dosage genotype file EUR_subset.dosage2.chr22last100.gz

Setting number of threads to 4
fam: EUR_subset.fam
bim(s): EUR_subset.bim
bed(s): EUR_subset.bed

=== Reading genotype data ===
[...snip]

Total elapsed time for analysis = 36.6498 sec

[user@cn0852 example]$ **exit**
exit
salloc.exe: Relinquishing job allocation 11198398
salloc.exe: Job allocation 11198398 has been revoked.

[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g., run\_example.sh in the example directory). For example:



```

#!/bin/bash
module load BOLT-LMM
bolt \
    --bfile=EUR_subset \
    --remove=EUR_subset.remove \
    --exclude=EUR_subset.exclude \
    --phenoFile=EUR_subset.pheno.covars \
    --phenoCol=PHENO \
    --covarFile=EUR_subset.pheno.covars \
    --covarCol=CAT_COV \
    --qCovarCol=QCOV{1:2} \
    --modelSnps=EUR_subset.modelSnps \
    --lmm \
    --LDscoresFile=../tables/LDSCORE.1000G_EUR.tab.gz \
    --numThreads=2 \
    --statsFile=example.stats \
    --dosageFile=EUR_subset.dosage.chr17first100 \
    --dosageFile=EUR_subset.dosage.chr22last100.gz \
    --dosageFidIidFile=EUR_subset.dosage.indivs \
    --statsFileDosageSnps=example.dosageSnps.stats \
    --impute2FileList=EUR_subset.impute2FileList.txt \
    --impute2FidIidFile=EUR_subset.impute2.indivs \
    --statsFileImpute2Snps=example.impute2Snps.stats \
    --dosage2FileList=EUR_subset.dosage2FileList.txt \
    --statsFileDosage2Snps=example.dosage2Snps.stats \
    2>&1 | tee example.log # log output written to stdout and stderr

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] run_example.sh
```

Note: For test runs: Copy /usr/local/apps/BOLT-LMM/v2.3/example and /usr/local/apps/BOLT-LMM/v2.3/tables to your local directory. Under the example directory, run\_example.sh and run\_example\_reml2.sh are batch scripts that can be used in slurm 


Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. BOLT-LMM.swarm). For example:



```

bolt --bfile=AFR_subset ...
bolt --bfile=AMR_subset ...
bolt --bfile=EAS_subset ...
bolt --bfile=EUR_subset ...
bolt --bfile=SAS_subset ...

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f BOLT-LMM.swarm [-g #] [-t #] --module BOLT-LMM
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g #  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t #  Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module BOLT-LMM  Loads the BOLT-LMM module for each subjob in the swarm 
 | |
 | |
 | |








