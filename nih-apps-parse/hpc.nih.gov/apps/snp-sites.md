

document.querySelector('title').textContent = 'SNP-SITES on Biowulf';
SNP-SITES on Biowulf


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



SNP-sites extracts single nucleotide polymorphisms (SNPs) from a multi-FASTA alignment outputs results in multiple formats for downstream analysis.



### References:


* [SNP-sites: rapid efficient extraction of SNPs from multi-FASTA alignments](https://doi.org/10.1099/mgen.0.000056). Andrew J. Page, Ben Taylor, Aidan J. Delaney, Jorge Soares, Torsten Seemann, Jacqueline A. Keane, Simon R. Harris, Microbial Genomics 2(4), 2016.


Documentation
* [SNP-sites Main Site](http://sanger-pathogens.github.io/snp-sites/)


Important Notes
* Module Name: snp-sites (see [the modules page](/apps/modules.html) for more information)



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load snp-sites**
[+] Loading singularity  3.4.2  on cn3103 
[+] Loading snp-sites 2.4.1  ... 

[user@cn3144 ~]$ **cd /data/user/SNP-SITES\_TEST**

[user@cn3144 ~]$ **snp-sites PMEN1.aln.gz > PMEN1.out**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. snp-sites.sh). For example:



```

#!/bin/bash
set -e
module load snp-sites
snp-sites /data/user/SNP-SITES_TEST/PMEN1.aln.gz > /data/user/SNP-SITES_TEST/PMEN1.out

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] TEMPLATE.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. snp-sites.swarm). For example:



```

snp-sites /data/user/SNP-SITES_TEST/PMEN1.aln.gz > /data/user/SNP-SITES_TEST/PMEN1.out
snp-sites /data/user/SNP-SITES_TEST/ST239.aln.gz > /data/user/SNP-SITES_TEST/ST239.out

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f snp-sites.swarm [-g #] [-t #] --module snp-sites
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module snp-sites Loads the snp-sites module for each subjob in the swarm 
 | |
 | |
 | |








