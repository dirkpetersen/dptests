

document.querySelector('title').textContent = "flanker";
Flanker: Gene-flank analysis tool


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



Flanker is a tool for studying the homology of gene-flanking sequences. It will annotate FASTA/multi-FASTA files for specified genes, then write the flanking sequences to new FASTA files. There is also an optional step to cluster the flanks by sequence identity.



### References:


* Blow J.
 [**Flanker: a tool for comparative genomics of gene flanking regions. Microbial Genomics. 2021**](https://doi.org/10.1099/mgen.0.000634)
*J Mol Biol. 2012 Jan 13;415(2):406-18.*


Documentation
* [Documentation](https://flanker.readthedocs.io/en/latest/)


Important Notes
* Module Name: flanker* Multithreaded/singlethreaded/MPI...
 * This application produces HTML reports. You can use [hpcdrive](/docs/hpcdrive.html) to view these reports on your local workstation.
* Environment variables set 
	+ TEMPLATE\_HOME* Example files in ???* Reference data in /fdb/TEMPLATE/



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

[user@cn3144 ~]$ **module load flanker**

[user@cn3144 ~]$  **flanker -h**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. TEMPLATE.sh). For example:



```

#!/bin/bash
set -e
module load flanker
flanker < TEMPLATE.in > TEMPLATE.out

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] TEMPLATE.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. TEMPLATE.swarm). For example:



```

TEMPLATE < TEMPLATE.in > TEMPLATE.out
TEMPLATE < TEMPLATE.in > TEMPLATE.out
TEMPLATE < TEMPLATE.in > TEMPLATE.out
TEMPLATE < TEMPLATE.in > TEMPLATE.out

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f TEMPLATE.swarm [-g #] [-t #] --module TEMPLATE
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module TEMPLATE Loads the TEMPLATE module for each subjob in the swarm 
 | |
 | |
 | |








