

document.querySelector('title').textContent = 'angsd on Biowulf';
angsd on Biowulf


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



 ANGSD is a software for analyzing next generation sequencing data. The software can handle a number of different input types from mapped reads to imputed genotype probabilities. Most methods take genotype uncertainty into account instead of basing the analysis on called genotypes. This is especially useful for low and medium depth data. The software is written in C++ and has been used on large sample sizes. 


Documentation
* [angsd Main Site](http://www.popgen.dk/angsd/index.php/ANGSD)


Important Notes
* Module Name: angsd (see [the modules page](/apps/modules.html) for more information)
* Singlethreaded app
* $ANGSD\_TESTDATA enviornment variable for example files



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

[user@cn3144 ~]$ **mkdir angsd && cd angsd**

[user@cn3144 ~]$ **module load angsd**

[user@cn3144 ~]$ **cp $ANGSD\_TESTDATA/\* .**

[user@cn3144 ~]$ **tar xf bams.tar.gz** 

[user@cn3144 ~]$ **for i in bams/\*.bam;do samtools index $i;done**

[user@cn3144 ~]$ **ls bams/\*.bam > bam.filelist**

[user@cn3144 ~]$ **angsd -b bam.filelist -GL 1 -doMajorMinor 1 -doMaf 2 -P 5**
  -> angsd version: 0.929 (htslib: 1.9-118-g2da4c7d-dirty) build(Feb 22 2019 11:59:23)
  -> No '-out' argument given, output files will be called 'angsdput'
[bammer_main] 10 samples in 10 input files
  -> Parsing 10 number of samples 
  -> Printing at chr: 20 pos:14078459 chunknumber 3400 contains 584 sites
  -> Done reading data waiting for calculations to finish
  -> Done waiting for threads
  -> Output filenames:
    ->"angsdput.arg"
    ->"angsdput.mafs.gz"
  -> Fri Feb 22 13:19:44 2019
  -> Arguments and parameters for all analysis are located in .arg file
  -> Total number of sites analyzed: 1702624
  -> Number of sites retained after filtering: 1702593 
  [ALL done] cpu-time used =  163.22 sec
  [ALL done] walltime used =  56.00 sec

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. angsd.sh). For example:



```

#!/bin/bash
set -e
module load angsd
angsd -b bam.filelist -GL 1 -doMajorMinor 1 -doMaf 2 -P 5

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --mem=10g angsd.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. angsd.swarm). For example:



```

angsd -b bam.filelist1 -out out1 -GL 1 -doMajorMinor 1 -doMaf 2 -P 5
angsd -b bam.filelist2 -out out2 -GL 1 -doMajorMinor 1 -doMaf 2 -P 5
angsd -b bam.filelist3 -out out3 -GL 1 -doMajorMinor 1 -doMaf 2 -P 5

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f angsd.swarm -g 10 --module angsd
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module angsd Loads the angsd module for each subjob in the swarm 
 | |
 | |
 | |








