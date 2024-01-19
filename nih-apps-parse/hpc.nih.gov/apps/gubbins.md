

document.querySelector('title').textContent = 'GUBBINS on Biowulf';
GUBBINS on Biowulf


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



Gubbins (Genealogies Unbiased By recomBinations In Nucleotide Sequences) is an algorithm that iteratively identifies loci containing elevated densities of base substitutions while concurrently constructing a phylogeny based on the putative point mutations outside of these regions.



### References:


* [Rapid phylogenetic analysis of large samples of recombinant bacterial whole genome sequences using Gubbins](https://doi.org/10.1093/nar/gku1196). Croucher N. J., Page A. J., Connor T. R., Delaney A. J., Keane J. A., Bentley S. D., Parkhill J., Harris S.R. Nucleic Acids Research, Volume 43, Issue 3, 18 February 2015, Page e15.


Documentation
* [Gubbins Main Site](http://sanger-pathogens.github.io/gubbins/)


Important Notes
* Module Name: gubbins (see [the modules page](/apps/modules.html) for more information)



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

[user@cn3144 ~]$ **module load gubbins**
[+] Loading singularity  3.4.2  on cn3144 
[+] Loading gubbins 2.3.4  ... 

[user@cn3144 ~]$ **cd /data/user/GUBBINS\_TEST**

[user@cn3144 ~]$ **ls**
PMEN1.aln.gz  ST239.aln.gz

[user@cn3144 ~]$ **gunzip PMEN1.aln.gz**

[user@cn3144 ~]$ **gunzip ST239.aln.gz**

[user@cn3144 ~]$ **ls**
PMEN1.aln  ST239.aln

[user@cn3144 ~]$ **gubbins PMEN1.aln**

[user@cn3144 ~]$ **ls PMEN1.\***
PMEN1.aln				PMEN1.filtered_polymorphic_sites.phylip  PMEN1.recombination_predictions.embl
PMEN1.aln.seq.joint.txt			PMEN1.final_tree.tre			 PMEN1.recombination_predictions.gff
PMEN1.branch_base_reconstruction.embl	PMEN1.node_labelled.final_tree.tre	 PMEN1.summary_of_snp_distribution.vcf
PMEN1.filtered_polymorphic_sites.fasta	PMEN1.per_branch_statistics.csv

[user@cn3144 ~]$ **gubbins ST239.aln**

[user@cn3144 ~]$ **ls ST239.\***
ST239.aln				ST239.filtered_polymorphic_sites.phylip  ST239.recombination_predictions.embl
ST239.aln.seq.joint.txt			ST239.final_tree.tre			 ST239.recombination_predictions.gff
ST239.branch_base_reconstruction.embl	ST239.node_labelled.final_tree.tre	 ST239.summary_of_snp_distribution.vcf
ST239.filtered_polymorphic_sites.fasta	ST239.per_branch_statistics.csv

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. gubbins.sh). For example:



```

#!/bin/bash
set -e
module load gubbins
cd /data/user/GUBBINS_TEST
gubbins PMEN1.aln

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] gubbins.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. gubbins.swarm). For example:



```

cd /data/user/GUBBINS_TEST; gubbins PMEN1.aln
cd /data/user/GUBBINS_TEST; gubbins ST239.aln

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f gubbins.swarm [-g #] [-t #] --module gubbins
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module gubbins Loads the gubbins module for each subjob in the swarm 
 | |
 | |
 | |








