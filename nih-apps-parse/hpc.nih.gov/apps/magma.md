

document.querySelector('title').textContent = 'MAGMA on Biowulf';
MAGMA on Biowulf


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



MAGMA is a tool for gene analysis and generalized gene-set analysis of GWAS data. It can be used to analyse both raw genotype data as well as summary SNP p-values from a previous GWAS or meta-analysis.



### References:


* [de Leeuw, Christiaan A., et al. "MAGMA: generalized gene-set analysis of GWAS data." *PLoS Comput Biol* 11.4 (2015): e1004219.](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004219)


Documentation
* [MAGMA Main Site](http://ctg.cncr.nl/software/magma)


Important Notes
* Module Name: MAGMA (see [the modules page](/apps/modules.html) for more information)
* Singlethreaded app
* Example files in /usr/local/apps/MAGMA/gene\_location/



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

[user@cn3144 ~]$ **module load MAGMA**

[user@cn3144 ~]$ **magma --annotate --snp-loc /usr/local/apps/MAGMA/gene\_location/example.snp.txt \**
**--gene-loc /usr/local/apps/MAGMA/gene\_location/NCBI37.3/NCBI37.3.gene.loc \**
**--out temp.test**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226

[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. MAGMA.sh). For example:



```

#!/bin/bash
module load MAGMA
magma --annotate --snp-loc /usr/local/apps/MAGMA/gene_location/example.snp.txt \
--gene-loc /usr/local/apps/MAGMA/gene_location/NCBI37.3/NCBI37.3.gene.loc \
--out temp.test

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch MAGMA.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. MAGMA.swarm). For example:



```

magma annotate --snp-loc snp1.txt --gene-loc NCBI37.3.gene.loc --out out1
magma annotate --snp-loc snp2.txt --gene-loc NCBI37.3.gene.loc --out out2
magma annotate --snp-loc snp3.txt --gene-loc NCBI37.3.gene.loc --out out3

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f MAGMA.swarm [-g #] [-t #] --module MAGMA
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module MAGMA Loads the MAGMA module for each subjob in the swarm 
 | |
 | |
 | |








