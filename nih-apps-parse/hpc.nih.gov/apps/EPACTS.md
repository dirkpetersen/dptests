

document.querySelector('title').textContent = 'EPACTS on Biowulf';
EPACTS on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
 |



EPACTS (Efficient and Parallelizable Association Container Toolbox) is a versatile software pipeline to perform various statistical tests for 
identifying genome-wide association from sequence data through a user-friendly interface, both to scientific analysts and to method developers.




### References:


* EPACTS was developed in the Abecasis lab at the University of Michigan. [[EPACTS webpage](http://genome.sph.umich.edu/wiki/EPACTS)]


Documentation
* [EPACTS Wiki](http://genome.sph.umich.edu/wiki/EPACTS)


Important Notes
* Module Name: EPACTS (see [the modules page](/apps/modules.html) for more information)
* Singlethreaded



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

[user@cn3144 ~]$ **module load EPACTS**

[user@cn3144 ~]$ **cp ${EPACTS\_DIR}/1000G\_\* .**

[user@cn3144 ~]$ **epacts single --vcf 1000G\_exome\_chr20\_example\_softFiltered.calls.vcf.gz \
 --ped 1000G\_dummy\_pheno.ped --min-maf 0.001 --chr 20 --pheno DISEASE \
 --cov AGE --cov SEX --test b.score --anno --out test --run 2**
Detected phenotypes with 2 unique values - 1 and 2 - considering them as binary phenotypes... re-encoding them into 1 and 2
Successfully written phenotypes and 2 covariates across 266 individuals
Processing chromosome 20...
Finished generating EPACTS Makefile
Running 2 parallel jobs of EPACTS
forkExecWait(): make -f /data/user/test.Makefile -j 2
Rscript /usr/local/share/EPACTS/epactsSingle.R --vanilla /usr/local /data/user/test.phe /data/user/test.cov /data/user/test.ind /data/user/1000G_exome_chr20_example_softFiltered.calls.vcf.gz 20:1-10000000 /data/user/test.20.1.10000000.epacts GT 0.001 1 3 1000000000 0.5 0 FALSE single.b.score
Rscript /usr/local/share/EPACTS/epactsSingle.R --vanilla /usr/local /data/user/test.phe /data/user/test.cov /data/user/test.ind /data/user/1000G_exome_chr20_example_softFiltered.calls.vcf.gz 20:10000001-20000000 /data/user/test.20.10000001.20000000.epacts GT 0.001 1 3 1000000000 0.5 0 FALSE single.b.score
Loading required package: epactsR
Loading required package: epactsR
NOTICE - Reading VCF took 1 seconds
[....]
zcat /data/user/test.epacts.gz | awk '$9 != "NA" { print $0 }' | sort -g -k 9 | head -n 5000 > /data/user/test.epacts.top5000
touch /data/user/test.epacts.OK

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. EPACTS.sh). For example:



```

#!/bin/bash
set -e

cd /data/$USER
module load EPACTS
cp ${EPACTS_DIR}/1000G_* .
epacts single --vcf  1000G_exome_chr20_example_softFiltered.calls.vcf.gz   \
       --ped  1000G_dummy_pheno.ped    --min-maf 0.001 --chr 20 --pheno DISEASE \
       --cov AGE --cov SEX --test b.score --anno  --out test --run 2

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--mem=#] EPACTS.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. EPACTS.swarm). For example:



```

epacts single --vcf input1.vcf --ped pheno.ped --out out1 --run 2
epacts single --vcf input2.vcf --ped pheno.ped --out out2 --run 2
epacts single --vcf input3.vcf --ped pheno.ped --out out3 --run 2
epacts single --vcf input4.vcf --ped pheno.ped --out out4 --run 2

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f EPACTS.swarm [-g #] [-t #] --module EPACTS
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module EPACTS Loads the EPACTS module for each subjob in the swarm 
 | |
 | |
 | |








