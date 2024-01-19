

document.querySelector('title').textContent = 'IonQuant on Biowulf';
IonQuant on Biowulf


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



IonQuant is a fast and comprehensive tool for MS1 precursor intensity-based quantification for timsTOF PASEF DDA and non-timsTOF (e.g., Orbitrap) data. It enables label-free quantification with false discovery (FDR) controlled match-between-runs (MBR). It can also be used for quantification in labelling-based experiments such as those involving SILAC, dimethyl, or similar labelling strategies. IonQuant is available as part of [FragPipe](https://hpc.nih.gov/apps/fragpipe.html) (recommended option), but can also be run as a command-line tool.


### Reference:


* [Yu, Fengchao., et al. "Fast quantitative analysis of timsTOF PASEF data with MSFragger and IonQuant." *Molecular & Cell Proteomics* 19 (2020): 1575-1585.](https://doi.org/10.1074/mcp.tir120.002048)


Documentation
* [ionquant Main Site](http://ionquant.nesvilab.org/)


Important Notes
* Module Name: ionquant (see [the modules page](/apps/modules.html) for more information)
 * Multithreaded application. Use --num\_threads=N to specify the number of threads to run.
 * Limit memory usage with the java heap size option -Xmx (e.g. -Xmx3700m for 3700 MB or -Xmx32g for 32 GB).
 * This application can be executed at the command line or invoked through the [fragpipe](fragpipe.html) GUI frontend.
* Environment variables set 
	+ IONQUANT\_HOME
	+ IONQUANT\_JAR



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

[user@cn0852 ~]$ **module load ionquant**
[+] Loading ionquant  1.8.0  on cn0852
[+] Loading java 12.0.1  ...

[user@cn0852 ~]$ **java -jar $IONQUANT\_JAR -h | head**
IonQuant version IonQuant-1.8.0
Batmass-IO version 1.25.5
timsdata library version timsdata-2-7-0
(c) University of Michigan
System OS: Linux, Architecture: amd64
Java Info: 17.0.3.1, Java HotSpot(TM) 64-Bit Server VM, Oracle Corporation
JVM started with 2 GB memory
Usage:
        java -jar IonQuant.jar  --specdir  --psm  --psm ...
 OR

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. ionquant.sh). For example:



```

#!/bin/bash
set -e
module load ionquanti
cd /path/to/data/and/conf

java -Xmx${SLURM_MEM_PER_NODE}m -jar $IONQUANT_JAR --threads=${SLURM_CPUS_PER_TASK} --filelist <files>


```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] ionquant.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. ionquant.swarm). For example:



```

java -Xmx${SLURM_MEM_PER_NODE}m -jar $IONQUANT_JAR --num_threads=${SLURM_CPUS_PER_TASK} --filelist <files 1>
java -Xmx${SLURM_MEM_PER_NODE}m -jar $IONQUANT_JAR --num_threads=${SLURM_CPUS_PER_TASK} --filelist <files 2>
java -Xmx${SLURM_MEM_PER_NODE}m -jar $IONQUANT_JAR --num_threads=${SLURM_CPUS_PER_TASK} --filelist <files 3>
java -Xmx${SLURM_MEM_PER_NODE}m -jar $IONQUANT_JAR --num_threads=${SLURM_CPUS_PER_TASK} --filelist <files 4>

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f ionquant.swarm [-g #] [-t #] --module ionquant
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








