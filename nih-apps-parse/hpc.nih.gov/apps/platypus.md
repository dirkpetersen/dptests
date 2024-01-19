

document.querySelector('title').textContent = 'Platypus on Biowulf';
Platypus on Biowulf


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



Platypus is a tool designed for efficient and accurate variant-detection in high-throughput sequencing data.



### References:


* Andy Rimmer, Hang Phan, Iain Mathieson, Zamin Iqbal, Stephen R. F. Twigg, WGS500 Consortium, Andrew O. M. Wilkie, Gil McVean, Gerton Lunter. Integrating mapping-, assembly- and haplotype-based approaches for calling variants in clinical sequencing applications. Nature Genetics (2014) [doi:10.1038/ng.3036](https://dx.doi.org/10.1038/ng.3036)


Documentation
* [Platypus main site](http://www.well.ox.ac.uk/platypus)


Important Notes
* Module Name: platypus (see [the modules page](/apps/modules.html) for more information)



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --gres lscratch:1**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job
[user@cn3144 ~]$ **cd /lscratch/$SLURM\_JOBID**
[user@cn3144 46116226]$ **module load platypus**
[+] Loading platypus, version 0.8.1...
[user@cn3144 46116226]$ **platypus callVariants \**
> **--bamFiles /fdb/app\_testdata/bam/hg19/gcat\_set\_053.bam \**
> **--refFile /fdb/igenomes/Homo\_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa \**
> **--output gcat\_set\_053.vcf** 
2017-09-29 14:53:42,435 - INFO - Beginning variant calling
2017-09-29 14:53:42,435 - INFO - Output will go to gcat_set_053.vcf
2017-09-29 14:53:42,710 - INFO - Processing region chr1:0-100000. (Only printing this message every 10 regions of size 100000)
2017-09-29 14:53:43,061 - INFO - Processing region chr1:1000000-1100000. (Only printing this message every 10 regions of size 100000)
...
...
2017-09-29 15:04:50,820 - INFO - Processing region chrY:58300000-58400000. (Only printing this message every 10 regions of size 100000)
2017-09-29 15:04:51,941 - INFO - Processing region chrY:59300000-59373566. (Only printing this message every 10 regions of size 100000)
2017-09-29 15:04:51,947 - INFO - Merging output VCF file(s) into final file gcat_set_053.vcf
2017-09-29 15:04:55,587 - INFO - Finished merging VCF file(s)
2017-09-29 15:04:55,591 - INFO - Finished variant calling
[user@cn3144 46116226]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. platypus.sh). For example:



```

#!/bin/sh

set -e
module load platypus

platypus callVariants \
 --bamFiles /fdb/app_testdata/bam/hg19/gcat_set_053.bam \
 --refFile /fdb/igenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa \
 --output gcat_set_053.vcf

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] platypus.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. platypus.swarm). For example:



```

platypus callVariants --bamFiles input1.bam --refFile hg19.fa --output input1.vcf
platypus callVariants --bamFiles input2.bam --refFile hg19.fa --output input2.vcf
platypus callVariants --bamFiles input3.bam --refFile hg19.fa --output input3.vcf
platypus callVariants --bamFiles input4.bam --refFile hg19.fa --output input4.vcf

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f platypus.swarm [-g #] [-t #] --module platypus
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g #  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t #  Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module platypus  Loads the platypus module for each subjob in the swarm 
 | |
 | |
 | |








