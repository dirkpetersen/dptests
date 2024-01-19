

document.querySelector('title').textContent = 'Shapeit on HPC';
Shapeit on HPC


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

  **SHAPEIT** is a fast and accurate method for [estimation 
 of haplotypes](http://en.wikipedia.org/wiki/Haplotype_estimation) (aka phasing) from genotype or sequencing data.  

  

**SHAPEIT** has several notable features:


* Linear complexity with the number of SNPs and conditioning haplotypes.
* Whole chromosome GWAS scale datasets can be phased in a single run.
* Phasing individuals with any level of relatedness
* Phasing is multi-threaded to tailor computational times to your resources
* Handles X chromosomes phasing
* Phasing using a reference panel (eg.1,000 Genomes) to aid phasing
* Ideal for pre-phasing imputation together with IMPUTE2


### References:


* **SHAPEIT** has primarily been developed by [Dr 
 Olivier Delaneau](http://funpopgen.unige.ch/members/olivier-delaneau) through a collaborative project between the research 
 groups of [Prof Jean-Francois 
 Zagury](http://gba.cnam.fr/index_eng.php?index=jfz)at [CNAM](http://www.cnam.fr/) and [Prof 
 Jonathan Marchini](http://www.stats.ox.ac.uk/~marchini/)at [Oxford](http://www.ox.ac.uk/). Funding 
 for this project has been received from several sources : [CNAM](http://www.cnam.fr/), 
 [Peptinov](http://www.peptinov.fr/), [MRC](http://www.mrc.ac.uk/index.htm), 
 [Leverhulme](http://www.leverhulme.ac.uk/), [The 
 Wellcome Trust](http://www.wellcome.ac.uk/).


Documentation
* <https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html>



Important Notes
* Module Name: shapeit (see [the 
 modules page](/apps/modules.html) for more information)
* Multithreaded
* Example files in /usr/local/apps/shapeit/example





Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --cpus-per-task=4**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load shapeit**[user@cn3144 ~]$ **shapeit --input-bed gwas.bed gwas.bim gwas.fam --input-map genetic\_map.txt --output-max gwas.phased.haps gwas.phased.sample --thread $SLURM\_CPUS\_PER\_TASK**

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
module load shapeit
shapeit --input-bed gwas.bed gwas.bim gwas.fam --input-map genetic_map.txt --output-max gwas.phased.haps gwas.phased.sample --thread $SLURM_CPUS_PER_TASK
```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=4 [--mem=#] batch.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. job.swarm). For example:



```

cd dir1; shapeit --input-bed gwas.bed gwas.bim gwas.fam --input-map genetic_map.txt --output-max gwas.phased.haps gwas.phased.sample --thread $SLURM_CPUS_PER_TASK
cd dir2; shapeit --input-bed gwas.bed gwas.bim gwas.fam --input-map genetic_map.txt --output-max gwas.phased.haps gwas.phased.sample --thread $SLURM_CPUS_PER_TASK
cd dir3; shapeit --input-bed gwas.bed gwas.bim gwas.fam --input-map genetic_map.txt --output-max gwas.phased.haps gwas.phased.sample --thread $SLURM_CPUS_PER_TASK

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f job.swarm [-g #] -t 4 --module shapeit
```

where
 

|  |  |
| --- | --- |
| -g *#*  | Number of Gigabytes of memory required for each process (1 line in the swarm command file)
  |
| -t *#* | Number of threads/CPUs required for each process (1 line in the swarm command file).
  |
| --module  | Loads the module for each subjob in the swarm  |




