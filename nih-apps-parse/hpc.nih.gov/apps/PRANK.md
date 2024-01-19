

document.querySelector('title').textContent = 'PRANK on Biowulf';
PRANK on Biowulf


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



PRANK is a probabilistic multiple alignment program for DNA, codon and amino-acid sequences. PRANK is based on a novel algorithm that treats insertions correctly and avoids over-estimation of the number of deletion events. In addition, PRANK borrows ideas from maximum likelihood methods used in phylogenetics and correctly takes into account the evolutionary distances between sequences. Lastly, PRANK allows for defining a potential structure for sequences to be aligned and then, simultaneously with the alignment, predicts the locations of structural units in the sequences.



### References:


* Veidenberg A, Medlar A, Löytynoja A.
 [**Wasabi: An Integrated Platform for Evolutionary Sequence Analysis and Data Visualization.**](https://www.ncbi.nlm.nih.gov/pubmed/26635364)
*Mol Biol Evol. 2016 Apr;33(4):1126-30.*


Documentation
* [PRANK Homepage](http://wasabiapp.org/software/prank/)* Type **man prank** after loading the module


Important Notes
* Module Name: PRANK (see [the modules page](/apps/modules.html) for more information)
* Singlethreaded
* Environment variables set 
	+ PRANK\_HOME
	+ PRANK\_EXAMPLES* Example files in $PRANK\_EXAMPLES



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

[user@cn3144 ~]$ **module load PRANK**
[user@cn3144 ~]$ **prank -d=input\_file -t=tree\_file -o=output\_file -F -showxml**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. PRANK.sh). For example:



```

#!/bin/bash
module load PRANK
prank -d=input_file -t=tree_file -o=output_file -F -showxml

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] PRANK.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. PRANK.swarm). For example:



```

prank -d=input_file_1 -t=tree_file -o=output_file_1 -F -showxml
prank -d=input_file_2 -t=tree_file -o=output_file_2 -F -showxml
prank -d=input_file_3 -t=tree_file -o=output_file_3 -F -showxml
prank -d=input_file_4 -t=tree_file -o=output_file_4 -F -showxml

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f PRANK.swarm [-g #] [-t #] --module PRANK
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module PRANK Loads the PRANK module for each subjob in the swarm 
 | |
 | |
 | |








