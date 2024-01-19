

document.querySelector('title').textContent = 'mixer on Biowulf';
mixer on Biowulf


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



MiXeR is Causal Mixture Model for GWAS summary statistics. The version(1.3) installed here contains a Python port of MiXeR, wrapping the same C/C++ core. 
Also data preprocessing code sumstats.py is included too.






### References:


* D. Holland et al.*Beyond SNP Heritability: Polygenicity and Discoverability Estimated for Multiple Phenotypes with a Univariate Gaussian Mixture Model* PLOS Genetics. 2020
 [PubMed](https://pubmed.ncbi.nlm.nih.gov/32427991/) | 
 [Journal](https://doi.org/10.1371/journal.pgen.1008612)


Documentation
* mixer Github:[Github](https://github.com/precimed/mixer)


Important Notes
* Module Name: mixer (see [the modules page](/apps/modules.html) for more information)
 * mixer is installed as a container on Biowulf, so please follow the instruction in this webpage regarding the path to mix.py, libbgmg.so et al.
 
```

	python /opt/mixer/precimed/mixer.py --help
	
```




Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive --cpus-per-task=10 --mem=10G**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144]$ **module load mixer**
[user@cn3144]$ **mkdir /data/$USER/mixer\_test/**
[user@cn3144]$ **cd /data/$USER/mixer\_test/**
[user@cn3144]$ **cp -r ${MIXER\_TEST\_DATA:-none}/\* .**
[user@cn3144]$ **python /opt/python\_convert/sumstats.py csv --sumstats GWAS\_EA\_excl23andMe.txt.gz --out SSGAC\_EDU\_2018\_no23andMe.csv --force --auto --head 5 --n-val 766345**

***********************************************************************
* sumstats.py: utilities for GWAS summary statistics
* Version 1.0.0
* (C) 2016-2018 Oleksandr Frei and Alexey A. Shadrin
* Norwegian Centre for Mental Disorders Research / University of Oslo
* GNU General Public License v3
***********************************************************************
Call:
./sumstats.py csv \
	--sumstats GWAS_EA_excl23andMe.txt.gz \
	--out SSGAC_EDU_2018_no23andMe.csv \
	--force  \
	--auto  \
	--head 5 \
	--n-val 766345.0

...

Analysis finished at Thu Jul 29 10:32:54 2021
Total time elapsed: 1.0m:57.59s

[user@cn3144]$ **python /opt/python\_convert/sumstats.py zscore --sumstats SSGAC\_EDU\_2018\_no23andMe.csv | \
> python /opt/python\_convert/sumstats.py qc --exclude-ranges 6:26000000-34000000 --out SSGAC\_EDU\_2018\_no23andMe\_noMHC.csv --force**
***********************************************************************
* sumstats.py: utilities for GWAS summary statistics
* Version 1.0.0
* (C) 2016-2018 Oleksandr Frei and Alexey A. Shadrin
* Norwegian Centre for Mental Disorders Research / University of Oslo
* GNU General Public License v3
***********************************************************************
Call:
./sumstats.py qc \
	--out SSGAC_EDU_2018_no23andMe_noMHC.csv \
	--force  \
	--exclude-ranges ['6:26000000-34000000']

...

Analysis finished at Thu Jul 29 10:38:04 2021
Total time elapsed: 3.0m:5.710000000000008s
[user@cn3144]$ **python /opt/mixer/precimed/mixer.py ld \
 --lib /opt/mixer/src/build/lib/libbgmg.so \
 --bfile 1000G\_EUR\_Phase3\_plink/1000G.EUR.QC.22 \
 --out 1000G\_EUR\_Phase3\_plink/1000G.EUR.QC.22.run4.ld \
 --r2min 0.05 --ldscore-r2min 0.05 --ld-window-kb 30000**

INFO:root:__init__(lib_name=/opt/mixer/src/build/lib/libbgmg.so, context_id=0)
INFO:root:init_log(1000G_EUR_Phase3_plink/1000G.EUR.QC.22.run4.ld.log)
INFO:root:log_message(***********************************************************************
* mixer.py: Univariate and Bivariate Causal Mixture for GWAS
* Version 1.2.0
* (c) 2016-2020 Oleksandr Frei, Alexey A. Shadrin, Dominic Holland
* Norwegian Centre for Mental Disorders Research / University of Oslo
* Center for Multimodal Imaging and Genetics / UCSD
* GNU General Public License v3
***********************************************************************
Call:
./mixer.py ld \
	--out 1000G_EUR_Phase3_plink/1000G.EUR.QC.22.run4.ld \
	--lib /opt/mixer/src/build/lib/libbgmg.so \
	--bfile 1000G_EUR_Phase3_plink/1000G.EUR.QC.22 \
	--ldscore-r2min 0.05 \
	--ld-window-kb 30000.0
)
INFO:root:__init__(lib_name=/opt/mixer/src/build/lib/libbgmg.so, context_id=0)
INFO:root:log_message(Done)

[user@cn3144]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. mixer.sh). For example:



```


#!/bin/bash

#SBATCH --job-name=mixer_run
#SBATCH --time=2:00:00

#SBATCH --partition=norm
#SBATCH --nodes=1
#SBATCH --array=1-20

cd /data/$USER/mixer_test
module load mixer
 
# The example here only use chr22, please change accordingly for your own dataset.

python /opt/mixer/precimed/mixer.py snps \
     --lib /opt/mixer/src/build/lib/libbgmg.so \
     --bim-file 1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim \
     --ld-file 1000G_EUR_Phase3_plink/1000G.EUR.QC.@.run4.ld \
     --out 1000G_EUR_Phase3_plink/1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.rep${SLURM_ARRAY_TASK_ID}.snps \
     --chr2use 22 \
     --maf 0.05 --subset 2000000 --r2 0.8 --seed 1


python /opt/mixer/precimed/mixer.py fit1 \
      --trait1-file SSGAC_EDU_2018_no23andMe_noMHC.csv.gz \
      --out SSGAC_EDU_2018_no23andMe_noMHC.fit.rep${SLURM_ARRAY_TASK_ID}\
      --extract 1000G_EUR_Phase3_plink/1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.rep${SLURM_ARRAY_TASK_ID}.snps \
      --bim-file 1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim \
      --ld-file 1000G_EUR_Phase3_plink/1000G.EUR.QC.@.run4.ld \
      --chr2use 22 \
      --lib  /opt/mixer/src/build/lib/libbgmg.so

python /opt/mixer/precimed/mixer.py test1 \
      --trait1-file SSGAC_EDU_2018_no23andMe_noMHC.csv.gz \
      --load-params-file SSGAC_EDU_2018_no23andMe_noMHC.fit.rep${SLURM_ARRAY_TASK_ID}.json \
      --out SSGAC_EDU_2018_no23andMe_noMHC.test.rep${SLURM_ARRAY_TASK_ID} \
      --bim-file 1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim \
      --ld-file 1000G_EUR_Phase3_plink/1000G.EUR.QC.@.run4.ld \
      --chr2use 22 \
      --lib  /opt/mixer/src/build/lib/libbgmg.so



```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=20 --mem=2g mixer.sh
```









