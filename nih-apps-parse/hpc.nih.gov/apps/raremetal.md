

document.querySelector('title').textContent = 'RAREMETAL on Biowulf';
RAREMETAL on Biowulf


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



RAREMETAL is a program that facilitates the meta-analysis of rare variants from genotype arrays or sequencing.



### Reference:


* [Feng, Shuang, et al. "RAREMETAL: fast and powerful meta-analysis for rare variants." *Bioinformatics* 30.19 (2014): 2828-2829.](https://academic.oup.com/bioinformatics/article/30/19/2828/2422168)


Documentation
* [RAREMETAL Main Site](https://genome.sph.umich.edu/wiki/RAREMETAL)
* [RAREMETAL GitHub repository](https://github.com/statgen/raremetal)


Important Notes
* Module Name: raremetal (see [the modules page](/apps/modules.html) for more information)
 * Environment variables set 
	+ RAREMETAL\_HOME



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf ~]$ **sinteractive --cpus-per-task=4 --mem=8G --gres=lscratch:10**
salloc.exe: Pending job allocation 42848391
salloc.exe: job 42848391 queued and waiting for resources
salloc.exe: job 42848391 has been allocated resources
salloc.exe: Granted job allocation 42848391
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3101 are ready for job
srun: error: x11: no local DISPLAY defined, skipping

[user@cn3101 ~]$ **cd /lscratch/$SLURM\_JOB\_ID**

[user@cn3101 42848391]$ **cp -r /fdb/raremetal/4.15.1/raremetal\_tutorial .**

[user@cn3101 42848391]$ **cd raremetal\_tutorial/**

[user@cn3101 raremetal_tutorial]$ **module load raremetal samtools**
[+] Loading raremetal  4.15.1  on cn3101
[+] Loading singularity  3.4.2  on cn3101
[+] Loading samtools 1.9  ...

[user@cn3101 raremetal_tutorial]$ **raremetalworker --ped example1.ped --dat example1.dat --vcf example1.vcf.gz \
 --traitName QT1 --inverseNormal --makeResiduals --kinSave --kinGeno \
 --prefix STUDY1**

RAREMETALWORKER 4.15.1 -- A Forerunner of RareMetal
          (c) 2012-2016 Shuang Feng, Dajiang Liu, Sai Chen, Goncalo Abecasis
[...snip...]

[user@cn3101 raremetal_tutorial]$ **raremetalworker --ped example2.ped --dat example2.dat --vcf example2.vcf.gz \
 --traitName QT1 --inverseNormal --makeResiduals --kinSave --kinGeno \
 --prefix STUDY2**

RAREMETALWORKER 4.15.1 -- A Forerunner of RareMetal
          (c) 2012-2016 Shuang Feng, Dajiang Liu, Sai Chen, Goncalo Abecass
[...snip...]

[user@cn3101 raremetal_tutorial]$ **bgzip STUDY1.QT1.singlevar.score.txt**

[user@cn3101 raremetal_tutorial]$ **tabix -c "#" -s 1 -b 2 -e 2 STUDY1.QT1.singlevar.score.txt.gz**

[user@cn3101 raremetal_tutorial]$ **bgzip STUDY1.QT1.singlevar.cov.txt**

[user@cn3101 raremetal_tutorial]$ **tabix -c "#" -s 1 -b 2 -e 2 STUDY1.QT1.singlevar.cov.txt.gz**

[user@cn3101 raremetal_tutorial]$ **bgzip STUDY2.QT1.singlevar.score.txt**

[user@cn3101 raremetal_tutorial]$ **tabix -c "#" -s 1 -b 2 -e 2 STUDY2.QT1.singlevar.score.txt.gz**

[user@cn3101 raremetal_tutorial]$ **bgzip STUDY2.QT1.singlevar.cov.txt**

[user@cn3101 raremetal_tutorial]$ **tabix -c "#" -s 1 -b 2 -e 2 STUDY2.QT1.singlevar.cov.txt.gz**

[user@cn3101 raremetal_tutorial]$ **cat >summaryfiles<<"EOF"
STUDY1.QT1.singlevar.score.txt.gz
STUDY2.QT1.singlevar.score.txt.gz
EOF**

[user@cn3101 raremetal_tutorial]$ **cat >covfiles<<"EOF"
STUDY1.QT1.singlevar.cov.txt.gz
STUDY2.QT1.singlevar.cov.txt.gz
EOF**

[user@cn3101 raremetal_tutorial]$ **raremetal --summaryFiles summaryfiles --covFiles covfiles \
 --groupFile group.file --SKAT --burden --MB --VT --longOutput \
 --tabulateHits --hitsCutoff 1e-05 --prefix COMBINED.QT1 --hwe 1.0e-05 \
 --callRate 0.95 -**

RAREMETAL 4.15.1 -- A Tool for Rare Variants Meta-Analyses for Quantitative Traits
          (c) 2012-2017 Shuang Feng, Dajiang Liu, Sai Chen, Goncalo Abecasis

[...snip...]

[user@cn3101 raremetal_tutorial]$ **ls**
COMBINED.QT1.meta.burden.results      example2.vcf.gz
COMBINED.QT1.meta.MB.results          example2.vcf.gz.tbi
COMBINED.QT1.meta.plots.pdf           group.file
COMBINED.QT1.meta.singlevar.results   STUDY1.Empirical.Kinship.gz
COMBINED.QT1.meta.SKAT_.results       STUDY1.QT1.additive.plots.pdf
COMBINED.QT1.meta.tophits.burden.tbl  STUDY1.QT1.singlevar.cov.txt.gz
COMBINED.QT1.meta.tophits.MB.tbl      STUDY1.QT1.singlevar.cov.txt.gz.tbi
COMBINED.QT1.meta.tophits.SKAT.tbl    STUDY1.QT1.singlevar.score.txt.gz
COMBINED.QT1.meta.tophits.VT.tbl      STUDY1.QT1.singlevar.score.txt.gz.tbi
COMBINED.QT1.meta.VT_.results         STUDY1.singlevar.log
COMBINED.QT1.raremetal.log            STUDY2.Empirical.Kinship.gz
command_to_use                        STUDY2.QT1.additive.plots.pdf
covfiles                              STUDY2.QT1.singlevar.cov.txt.gz
example1.dat                          STUDY2.QT1.singlevar.cov.txt.gz.tbi
example1.ped                          STUDY2.QT1.singlevar.score.txt.gz
example1.vcf.gz                       STUDY2.QT1.singlevar.score.txt.gz.tbi
example1.vcf.gz.tbi                   STUDY2.singlevar.log
example2.dat                          summaryfiles
example2.ped

[user@cn3101 raremetal_tutorial]$ **exit**
exit
salloc.exe: Relinquishing job allocation 42848391

[user@biowulf ~]$


```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. raremetal.sh). For example:



```

#!/bin/bash
set -e
module load raremetal
raremetal --summaryFiles summaryfiles --covFiles covfiles \
    --groupFile group.file --SKAT --burden --MB --VT --longOutput \
    --tabulateHits --hitsCutoff 1e-05 --prefix COMBINED.QT1 --hwe 1.0e-05 \
    --callRate 0.95 -

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] raremetal.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. raremetal.swarm). For example:



```

raremetal --summaryFiles summaryfiles --covFiles covfiles \
    --groupFile group.file --SKAT --burden --MB --VT --longOutput \
    --tabulateHits --hitsCutoff 1e-05 --prefix COMBINED.QT1 --hwe 1.0e-05 \
    --callRate 0.95 -
raremetal --summaryFiles summaryfiles --covFiles covfiles \
    --groupFile group.file --SKAT --burden --MB --VT --longOutput \
    --tabulateHits --hitsCutoff 1e-05 --prefix COMBINED.QT2 --hwe 1.0e-05 \
    --callRate 0.95 -
raremetal --summaryFiles summaryfiles --covFiles covfiles \
    --groupFile group.file --SKAT --burden --MB --VT --longOutput \
    --tabulateHits --hitsCutoff 1e-05 --prefix COMBINED.QT3 --hwe 1.0e-05 \
    --callRate 0.95 -

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f raremetal.swarm [-g #] [-t #] --module raremetal
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module raremetal Loads the RAREMETAL module for each subjob in the swarm 
 | |
 | |
 | |








