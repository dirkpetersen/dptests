

document.querySelector('title').textContent = 'Primer3 on Biowulf';
Primer3 on Biowulf


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



Primer3 is a program for designing PCR primers.
Primer3 can also design hybridization probes and sequencing primers.



### References:


* Untergasser A, Cutcutache I, Koressaar T, Ye J, Faircloth BC, Remm M and Rozen SG.
[**Primer3--new capabilities and interfaces.**](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3424584/)
*Nucleic Acids Res. 2012 Aug 1;40(15):e115*
* Koressaar T and Remm M.
 [**Enhancements and modifications of primer design program Primer3.**](https://www.ncbi.nlm.nih.gov/pubmed/17379693)
*Bioinformatics 2007;23(10):1289-1291.*
* Koressaar T, Lepamets M, Kaplinski L, Raime K, Andreson R and Remm M.
 [**Primer3\_masker: integrating masking of template sequence with primer design software.**](https://www.ncbi.nlm.nih.gov/pubmed/29360956)
*Bioinformatics 2018;34(11):1937-1938.*


Documentation
* [Primer3 Main Site](https://primer3.org/)


Important Notes
* Module Name: primer3 (see [the modules page](/apps/modules.html) for more information)
 * This application produces HTML reports. You can use [hpcdrive](/docs/hpcdrive.html) to view these reports on your local workstation.



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

[user@cn3144 ~]$ **module load primer3**

```


Following the [Primer3 manual's example](https://primer3.org/manual.html#example), create the file example with the following contents:

```

SEQUENCE_ID=example
SEQUENCE_TEMPLATE=GTAGTCAGTAGACNATGACNACTGACGATGCAGACNACACACACACACACAGCACACAGGTATTAGTGGGCCATTCGATCCCGACCCAAATCGATAGCTACGATGACG
SEQUENCE_TARGET=37,21
PRIMER_TASK=generic
PRIMER_PICK_LEFT_PRIMER=1
PRIMER_PICK_INTERNAL_OLIGO=1
PRIMER_PICK_RIGHT_PRIMER=1
PRIMER_OPT_SIZE=18
PRIMER_MIN_SIZE=15
PRIMER_MAX_SIZE=21
PRIMER_MAX_NS_ACCEPTED=1
PRIMER_PRODUCT_SIZE_RANGE=75-100
P3_FILE_FLAG=1
SEQUENCE_INTERNAL_EXCLUDED_REGION=37,21
PRIMER_EXPLAIN_FLAG=1
=

```


Then run primer3 and show the output files:


```

[user@cn3144 ~]$ **primer3\_core < example** 

[user@cn3144 ~]$ **ls**
example
example.for
example.int
example.rev

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. primer3.sh). For example:



```

#!/bin/bash
set -e
module load primer3
primer3_core < example

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] primer3.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. primer3.swarm). For example:



```

primer3_core < seq1
primer3_core < seq2
primer3_core < seq3

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f primer3.swarm [-g #] [-t #] --module primer3
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module primer3 Loads the primer3 module for each subjob in the swarm 
 | |
 | |
 | |












