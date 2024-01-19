

document.querySelector('title').textContent = 'msfragger on Biowulf';
msfragger on Biowulf


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



MSFragger is an ultrafast database search tool for peptide identification in mass spectrometry-based proteomics. It has demonstrated excellent performance across a wide range of datasets and applications. MSFragger is suitable for standard shotgun proteomics analyses as well as large datasets (including timsTOF PASEF data), enzyme unconstrained searches (e.g. peptidome), ‘open’ database searches (i.e. precursor mass tolerance set to hundreds of Daltons) for identification of modified peptides, and glycopeptide identification (N-linked and O-linked) with MSFragger Glyco mode.



### Reference:


* [Kong, Andy T., et al. "MSFragger: ultrafast and comprehensive peptide identification in mass spectrometry–based proteomics." *Nature methods* 14.5 (2017): 513-520.](https://www.nature.com/articles/nmeth.4256)


Documentation
* [msfragger Main Site](http://msfragger.nesvilab.org/)


Important Notes
* Module Name: msfragger (see [the modules page](/apps/modules.html) for more information)
 * Multithreaded application. Use --num\_threads=N to specify the number of threads to run.
 * Limit memory usage with the java heap size option -Xmx (e.g. -Xmx3700m for 3700 MB or -Xmx32g for 32 GB).
 * This application can be executed at the command line or invoked through the [fragpipe](fragpipe.html) GUI frontend.
* Environment variables set 
	+ MSFRAGGER\_HOME
	+ MSFRAGGER\_JAR



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf ~]$ **sinteractive --cpus-per-task=8 --mem=16g**
salloc.exe: Pending job allocation 4852341
salloc.exe: job 4852341 queued and waiting for resources
salloc.exe: job 4852341 has been allocated resources
salloc.exe: Granted job allocation 4852341
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn0852 are ready for job
srun: error: x11: no local DISPLAY defined, skipping

[user@cn0852 ~]$ **module load msfragger**
[+] Loading msfragger  3.1.1  on cn0852
[+] Loading java 12.0.1  ...

[user@cn0852 ~]$ **java -jar $MSFRAGGER\_JAR -h | head**
Usage:
        To perform a search either:
                1) java -jar MSFragger.jar  
 To generate default parameter files use --config flag. E.g. "java -jar MSFragger.jar --config"
 Or:
 2) java -jar MSFragger.jar  
Options:
--num\_threads  # Number of CPU threads to use, should be set to the
 # number of logical processors; A value of 0
 # (auto-detect) will cause MSFragger to use the

[user@cn0852 ~]$ **java -Xmx16g -jar $MSFRAGGER\_JAR --num\_threads=8 params.conf raw.d**
MSFragger version MSFragger-3.1.1
Batmass-IO version 1.19.5
timsdata library version timsdata-2-7-0
(c) University of Michigan
RawFileReader reading tool. Copyright (c) 2016 by Thermo Fisher Scientific, Inc. All rights reserved.
System OS: Linux, Architecture: amd64
Java Info: 12.0.1, Java HotSpot(TM) 64-Bit Server VM, Oracle Corporation
JVM started with 16 GB memory
[snip...]

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. msfragger.sh). For example:



```

#!/bin/bash
set -e
module load msfragger
cd /path/to/data/and/conf

java -Xmx${SLURM_MEM_PER_NODE}m -jar $MSFRAGGER_JAR --num_threads=${SLURM_CPUS_PER_TASK} params.conf raw.d


```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] msfragger.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. msfragger.swarm). For example:



```

java -Xmx${SLURM_MEM_PER_NODE}m -jar $MSFRAGGER_JAR --num_threads=${SLURM_CPUS_PER_TASK} params.conf raw1.d
java -Xmx${SLURM_MEM_PER_NODE}m -jar $MSFRAGGER_JAR --num_threads=${SLURM_CPUS_PER_TASK} params.conf raw2.d
java -Xmx${SLURM_MEM_PER_NODE}m -jar $MSFRAGGER_JAR --num_threads=${SLURM_CPUS_PER_TASK} params.conf raw3.d
java -Xmx${SLURM_MEM_PER_NODE}m -jar $MSFRAGGER_JAR --num_threads=${SLURM_CPUS_PER_TASK} params.conf raw4.d

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f msfragger.swarm [-g #] [-t #] --module msfragger
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








