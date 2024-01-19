

document.querySelector('title').textContent = 'Beagle on Biowulf';
Beagle on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
 |


BEAGLE is a software program for imputing genotypes, inferring haplotype phase, and performing genetic association analysis. BEAGLE is designed to analyze large-scale data sets with hundreds of thousands of markers genotyped on thousands of samples. BEAGLE can


* phase genotype data (i.e. infer haplotypes) for unrelated individuals, parent-offspring pairs, and parent-offspring trios.
* infer sporadic missing genotype data.
* impute ungenotyped markers that have been genotyped in a reference panel.
* perform single marker and haplotypic association analysis.


Beagle was developed by Brian Browning at the University of Washington.


Documentation
* Beagle Main Site: [Beagle at Univ. of Washington](http://faculty.washington.edu/browning/beagle/beagle.html)


Important Notes
* Module Name: Beagle (see [the modules page](/apps/modules.html) for more information)
 * Multithreaded
 * Unusual environment variables set:
	+ **$BEAGLEDIR** Beagle installation directory
	+ **$BEAGLE\_JAR** path to the Beagle jar file* Example files in $BEAGLEDIR



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --cpus-per-task=8 --mem=8g**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load Beagle**
[user@cn3144 ~]$ **java -jar $BEAGLE\_JAR nthreads=$SLURM\_CPUS\_PER\_TASK gt=test.08Jun17.d8b.vcf.gz out=out.gt**
beagle.08Jun17.d8b.jar (version 4.1)
Copyright (C) 2014-2015 Brian L. Browning
Enter "java -jar beagle.08Jun17.d8b.jar" for a summary of command line arguments.
Start time: 10:15 AM EDT on 31 Aug 2017

Command line: java -Xmx27305m -jar beagle.jar
  nthreads=8
  gt=test.08Jun17.d8b.vcf.gz
  out=out.gt

No genetic map is specified: using 1 cM = 1 Mb

reference samples:       0
target samples:        191

[ ... ]

Number of markers:                1356
Total time for building model: 6 seconds
Total time for sampling:       8 seconds
Total run time:                15 seconds

End time: 10:15 AM EDT on 31 Aug 2017
beagle.08Jun17.d8b.jar (version 4.1) finished
user@cn3144 ~]$
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. Beagle.sh). For example:



```

#!/bin/bash
module load Beagle
java -jar $BEAGLE_JAR nthreads=$SLURM_CPUS_PER_TASK gt=test.08Jun17.d8b.vcf.gz out=out.gt

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] Beagle.sh
```







