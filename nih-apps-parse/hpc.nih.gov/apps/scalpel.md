

document.querySelector('title').textContent = 'Scalpel on HPC';
Scalpel on HPC


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

  Scalpel is a software package for detecting INDELs (INsertions and DELetions) 
 mutations in a reference genome which has been sequenced with next-generation 
 sequencing technology (e.g., [Illumina](http://www.illumina.com/)). 
 Scalpel is designed to perform localized micro-assembly of specific regions 
 of interest with the goal of detecting mutations with high accuracy and 
 increased power. It is based on the de Bruijn graph assembly paradigm and 
 implements an on-the-fly repeat composition analysis coupled with a self-tuning 
 k-mer strategy to increase specificity in regions characterized by complex 
 repeat structures. It supports three different modes of operation:


* Single: in single mode scalpel detects INDELs in one single dataset 
 (e.g., one individual exome).
* Denovo: in denovo mode scalpel detects de novo INDELs in one family 
 of four individuals (mom, dad, aff, sib).
* Somatic: in somatic mode scalpel detects somatic INDELs for a tumor/sample 
 pair.


 For all the modes of operation, scalpel requires that the raw reads have 
 been previously aligned with BWA using default parameters. See BWA description 
 for more info.


### 


Documentation * <https://github.com/mjafin/scalpel>



Important Notes * Module Name: scalpel (see [the 
 modules page](/apps/modules.html) for more information)





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

[user@cn3144 ~]$ **module load scalpel**
[user@cn3144 ~]$ **scalpel-discovery --single --bam file.bam --bed regions.bed --ref genome.fa**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```




Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. batch.sh). For example:



```

#!/bin/bash
set -e
module load scalpel
scalpel-discovery --single --bam file.bam --bed regions.bed --ref genome.fa
```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--mem=#] batch.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. job.swarm). For example:



```

cd dir1; scalpel-discovery --single --bam file.bam --bed regions.bed --ref genome.fa
cd dir2; scalpel-discovery --single --bam file.bam --bed regions.bed --ref genome.fa
cd dir3; scalpel-discovery --single --bam file.bam --bed regions.bed --ref genome.fa

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f job.swarm [-g #] --module scalpel
```

where
 

|  |  |
| --- | --- |
| -g *#*  | Number of Gigabytes of memory required for each process 
 (1 line in the swarm command file)  |
| --module  | Loads the module for each subjob in the swarm  |




