

document.querySelector('title').textContent = "yak";
yak on Biowulf


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



Yak is initially developed for two specific use cases: 1) to robustly estimate the base accuracy of CCS reads and assembly contigs, and 2) to investigate the systematic error rate of CCS reads. It achieves the goals by comparing sequences to the k-mer spectrum of short reads or by comparing spectra. No reference genome or truth data is needed.




It is worth noting that estimating base accuracy is tricky. When the accuracy approaches Q50, both unsampled and erroneous k-mers in short reads may interfere with a naive estimator. Yak introduces an empirical model to address this issue. Its estimate is less affected by the coverage and the quality of short reads.




Documentation
* [yak on GitHub](https://github.com/lh3/yak)


Important Notes
* Module Name: yak (see [the modules page](/apps/modules.html) for more information)
 * Multithreaded. Use the -t flag to set the number of threads
* Environment variables set 
	+ YAK\_HOME



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive --mem 32g**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load yak**

[user@cn3144 ~]$ **yak count -b37 -o ccs.yak /fdb/app\_testdata/fastq/Homo\_sapiens/HG002\_NA24385\_son-PacBio\_CCS\_15kb/m54238\_180901\_011437.Q20.fastq.gz**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. yak.sh). For example:



```

#!/bin/bash
set -e
module load yak

yak count -b37 -o ccs.yak /fdb/app_testdata/fastq/Homo_sapiens/HG002_NA24385_son-PacBio_CCS_15kb/m54238_180901_011437.Q20.fastq.gz

yak inspect ccs.yak > ccs.hist

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] yak.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. yak.swarm). For example:



```

yak count -b37 -o sample1.yak sample1.fq.gz
yak count -b37 -o sample2.yak sample2.fq.gz
yak count -b37 -o sample3.yak sample3.fq.gz
yak count -b37 -o sample4.yak sample4.fq.gz

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f yak.swarm [-g #] [-t #] --module yak
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module yak Loads the yak module for each subjob in the swarm 
 | |
 | |
 | |








