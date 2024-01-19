

document.querySelector('title').textContent = 'MToolBox on Biowulf';
MToolBox on Biowulf


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



MToolBox is a highly automated bioinformatics pipeline to reconstruct and analyze human mitochondrial DNA from high throughput sequencing data. It includes an updated computational strategy to assemble mitochondrial genomes from whole exome and/or genome sequencing and an improved fragment-classify tool for haplogroup assignment, functional and prioritization analysis of mitochondrial variants. It also provides pathogenicity scores, profiles of genome variability and disease-associations for mitochondrial variants and a Variant Call Format file featuring allele-specific heteroplasmy.



### References:


* Calabrese C, Simone D, Diroma MA, Santorsola M, Gutt C, Gasparre G, Picardi E, Pesole G, Attimonelli M. MToolBox: a highly automated pipeline for heteroplasmy annotation and prioritization analysis of human mitochondrial variants in high-throughput sequencing. Bioinformatics. 2014 Jul 14. PMID: 25028726


Documentation
* [MToolBox Web Site](https://github.com/mitoNGS/MToolBox)


Important Notes
* Module Name: mtoolbox (see [the modules page](/apps/modules.html) for more information)



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

[user@cn3144 ~]$ **module load mtoolbox**
[user@cn3144 ~]$ **cd /path/to/mtoolbox/input/dir**
[user@cn3144 ~]$ **MToolbox.sh -i *input file*** 
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. mtoolbox.sh). For example:



```

#!/bin/bash
set -e
module load mtoolbox
cd *MToolbox/input/files/dir*
MToolBox.sh -i fastq -I -M -r RCRS

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] mtoolbox.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. mtoolbox.swarm). For example:



```

cd /MToolbox/input/dir1 && MToolbox.sh -i <input_format> -r <reference_sequence> -m "<mapExome_options>" -a "<assembleMTgenome_options>" -c "<mt-classifier_options>"
cd /MToolbox/input/dir2 && Mtoolbox.sh -i <input_format> -r <reference_sequence> -m "<mapExome_options>" -a "<assembleMTgenome_options>" -c "<mt-classifier_options>"
cd /MToolbox/input/dir3 && MToolbox.sh -i <input_format> -r <reference_sequence> -m "<mapExome_options>" -a "<assembleMTgenome_options>" -c "<mt-classifier_options>"
cd /MToolbox/input/dir4 && MToolbox.sh -i <input_format> -r <reference_sequence> -m "<mapExome_options>" -a "<assembleMTgenome_options>" -c "<mt-classifier_options>"

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f mtoolbox.swarm [-g #] [-t #] --module mtoolbox
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module mtoolbox Loads the MToolBox module for each subjob in the swarm 
 | |
 | |
 | |








