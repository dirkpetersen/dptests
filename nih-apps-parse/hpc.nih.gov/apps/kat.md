

document.querySelector('title').textContent = 'KAT on Biowulf';
KAT on Biowulf


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



 KAT (K-mer Analysis Toolkit) provides a suite of tools that, through the use of k-mer counts, help the user address or identify issues such as determining sequencing completeness for assembly, assessing sequencing bias, identifying contaminants, validating genomic assemblies and filtering content. KAT is geared primarily to work with high-coverage genomic reads from Illumina devices, although can work with any fasta or fastq sequence file.



### References:


* Daniel Mapleson, Gonzalo Garcia Accinelli, George Kettleborough, Jonathan Wright, and Bernardo J. Clavijo.
 KAT: A K-mer Analysis Toolkit to quality control NGS datasets and genome assemblies.
 Bioinformatics, 2016. doi: 
 [**10.1093/bioinformatics/btw663**](https://doi.org/10.1093/bioinformatics/btw663)


Documentation
* [KAT Main Site](https://www.earlham.ac.uk/kat-tools)
* [Documentation](https://kat.readthedocs.io/en/latest/)
* [KAT on GitHub](https://github.com/TGAC/KAT)


Important Notes

* Module Name: kat (see [the modules page](/apps/modules.html) for more information)




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

[user@cn3144 ~]$ **module load kat**

[user@cn3144 ~]$ **kat hist /fdb/app\_testdata/fastq/H\_sapiens/hg100\_1m\_pe?.fq.gz** 

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. KAT.sh). For example:



```

#!/bin/bash
set -e
module load kat
kat hist /fdb/app_testdata/fastq/H_sapiens/hg100_1m_pe?.fq.gz

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] KAT.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. kat.swarm). For example:



```

kat hist sample1.fastq
kat hist sample2.fastq
kat hist sample3.fastq
kat hist sample4.fastq

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f kat.swarm [-g #] [-t #] --module kat
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module kat Loads the KAT module for each subjob in the swarm 
 | |
 | |
 | |








