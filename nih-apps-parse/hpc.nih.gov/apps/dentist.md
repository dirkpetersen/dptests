

document.querySelector('title').textContent = 'DENTIST on HPC';
DENTIST on HPC


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

 DENTIST (Detecting Errors iN analyses of summary 
staTISTics) is a quality control (QC) tool for summary-level data from 
genome-wide association studies (GWASs). It leverages the difference between
the observed GWAS test-statistic of a variant and its predicted value (using
the neighbouring variants and linkage equilibrium (LD) data from a reference
panel) to remove problematic variants. It can detect genotyping/imputation
errors in either the original GWAS or the LD reference samples, allelic errors
(i.e., the effect alleles of the variants are mislabelled) in the GWAS summary
data, as well as heterogeneity between the GWAS and LD reference samples.

 ### 


References * [Chen, W., Wu, Y., Zheng, Z. *et al.* Improved analyses of GWAS summary statistics by reducing data heterogeneity and errors. *Nat Commun* **12**, 7117 (2021).](https://www.nature.com/articles/s41467-021-27438-7)


Documentation * <https://github.com/Yves-CHEN/DENTIST>



Important Notes * Module Name: dentist (see [the modules 
 page](/apps/modules.html) for more information)
* DENTIST can be run in a multi-threaded mode by using the
--thread-num <<>num<>>
 option. Only do this if you have
allocated multiple processors in your sinteractive call or your job script.
 




Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf ~]$ **sinteractive**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load dentist**
[user@cn3144 ~]$ **cd /data/$USER/dir**
[user@cn3144 dir]$ **dentist**
*******************************************************************
* DENTIST (Detecting Errors iN analyses of summary staTISTics)
* Version 1.1.0.0
* (C) 2018 Wenhan Chen, Zhihong Zhu and Jian Yang
* The University of Queensland
* MIT License
*******************************************************************
Flags include:
--gwas-summary,--bfile,--bld,--out,
--chrID,--thread-num,--dup-threshold,--p-value-threshold,
--GWAS-pvalue-threshold,--delta-MAF,--maf,--extract,
--target,--target-bp,--radius,--with-NA-geno,
--wind-dist,--wind,--debug,--load-LD,
--iteration-num,--LD-unit-in-byte,--SVD-trunc-prop,--check-LD,
--write-LD,--freq,--impute,
[user@cn3144 dir]$ **dentist --gwas-summary summary\_data --bfile ref --out prefix**

[user@cn3144 dir]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```




Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. batch.sh). For example:



```

#!/bin/bash
set -e
module load dentist 
dentist --gwas-summary summary_data --bfile ref --out prefix
```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--mem=#] batch.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. job.swarm). For example:



```

cd dir1; dentist --gwas-summary summary_data --bfile ref --out prefix
cd dir2; dentist --gwas-summary summary_data --bfile ref --out prefix
cd dir3; dentist --gwas-summary summary_data --bfile ref --out prefix

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f job.swarm [-g #] --module dentist
```

where
 

|  |  |
| --- | --- |
| -g *#*  | Number of Gigabytes of memory required for each process 
 (1 line in the swarm command file)  |
| --module  | Loads the module for each subjob in the swarm  |




