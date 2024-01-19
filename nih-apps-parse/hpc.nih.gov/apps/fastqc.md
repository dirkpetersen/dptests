

document.querySelector('title').textContent = 'FastQC on Biowulf';
FastQC on Biowulf


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

  FastQC aims to provide a simple way to do some quality control checks on raw 
 sequence data coming from high throughput sequencing pipelines. It provides 
 a modular set of analyses which you can use to give a quick impression of 
 whether your data has any problems of which you should be aware before doing 
 any further analysis.


The main functions of FastQC are


- Import of data from BAM, SAM or FastQ files (any variant)  

 - Providing a quick overview to tell you in which areas there may be problems  

 - Summary graphs and tables to quickly assess your data  

 - Export of results to an HTML based permanent report  

 - Offline operation to allow automated generation of reports without running 
 the interactive application  




### References:

 * FastQC is developed by [Simon 
 Andrews, Babraham Bioinformatics](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/).


Documentation * <http://www.bioinformatics.babraham.ac.uk/projects/fastqc/>


Important Notes * Module Name: fastqc (see [the modules 
 page](/apps/modules.html) for more information)



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

[user@cn3144 ~]$ **module load fastqc**
[user@cn3144 ~]$ **fastqc -o output\_dir [-f fastq|bam|sam] -c contaminant\_file seqfile1 .. seqfileN**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit). Create a batch input file (e.g. fastqc.sh). For example:



```

#!/bin/bash
set -e
module load fastqc
fastqc -o output_dir [-f fastq|bam|sam] -c contaminant_file seqfile1 .. seqfileN
```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --mem=10g fastqc.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. fastqc.swarm). For example:



```
cd dir1;fastqc -o output_dir [-f fastq|bam|sam] -c contaminant_file seqfile1 .. seqfileN
cd dir2;fastqc -o output_dir [-f fastq|bam|sam] -c contaminant_file seqfile1 .. seqfileN
cd dir3;fastqc -o output_dir [-f fastq|bam|sam] -c contaminant_file seqfile1 .. seqfileN
cd dir4;fastqc -o output_dir [-f fastq|bam|sam] -c contaminant_file seqfile1 .. seqfileN
cd dir5;fastqc -o output_dir [-f fastq|bam|sam] -c contaminant_file seqfile1 .. seqfileN

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f fastqc.swarm -g 10 --module fastqc
```

where
 

|  |  |  |  |
| --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in 
 the swarm command file) 
 | --module fastqc Loads the fastqc module for each subjob in the swarm 
  | |
 | |








