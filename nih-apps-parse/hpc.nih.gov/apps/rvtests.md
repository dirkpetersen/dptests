

document.querySelector('title').textContent = 'Rvtests on Biowulf';
Rvtests on Biowulf


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



Rvtests, which stands for Rare Variant tests, is a flexible software package for genetic association analysis for sequence datasets. Since its inception, rvtests was developed as a comprehensive tool to support genetic association analysis and meta-analysis. It can analyze both unrelated individual and related (family-based) individuals for both quantitative and binary outcomes. It includes a variety of association tests (e.g. single variant score test, burden test, variable threshold test, SKAT test, fast linear mixed model score test). It takes VCF/BGEN/PLINK format as genotype input file and takes PLINK format phenotype file and covariate file.


### References:


* Xiaowei Zhan, et al. RVTESTS: An Efficient and Comprehensive Tool for Rare Variant Association Analysis Using Sequence Data. Bioinformatics 2016 32: 1423-1426. doi:10.1093/bioinformatics/btw079


Documentation
* [Rvtests Main Site](https://github.com/zhanxw/rvtests)


Important Notes
* Module Name: rvtests (see [the modules page](/apps/modules.html) for more information)
* singlethreaded
* Example files in /usr/local/apps/rvtests/2.0.6/example/



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

[user@cn3144 ~]$ **module load rvtests**

[user@cn3144 ~]$ **RV\_DIR=/usr/local/apps/rvtests/2.0.6/example/**

[user@cn3144 ~]$ **rvtest --inVcf $RV\_DIR/example.vcf --pheno $RV\_DIR/pheno --out output --single wald,score**
Thank you for using rvtests (version: 20171009, git: Unknown)
  For documentations, refer to http://zhanxw.github.io/rvtests/
  For questions and comments, send to Xiaowei Zhan 
 For bugs and feature requests, please submit at: https://github.com/zhanxw/rvtests/issues
[...]
[INFO] Analysis begins with [ 9 ] samples...
[INFO] Impute missing genotype to mean (by default)
[INFO] Analysis started
[INFO] Analyzed [ 3 ] variants
[INFO] Analysis ends at: Mon Apr 9 11:50:48 2018
[INFO] Analysis took 2 seconds
RVTESTS finished successfully

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. rvtests.sh). For example:



```

#!/bin/bash
set -e
module load rvtests
RV_DIR=/usr/local/apps/rvtests/2.0.6/example/
rvtest --inVcf $RV_DIR/example.vcf --pheno $RV_DIR/pheno --out output --single wald,score
```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch rvtests.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. rvtests.swarm). For example:



```

rvtest --inVcf example.1vcf --pheno pheno1 --out output1 --single wald,score
rvtest --inVcf example.2vcf --pheno pheno2 --out output2 --single wald,score
rvtest --inVcf example.3vcf --pheno pheno3 --out output3 --single wald,score

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f rvtests.swarm --module rvtests
```

where


|  |  |
| --- | --- |
| --module rvtests Loads the rvtests module for each subjob in the swarm 
 | |








