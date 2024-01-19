

document.querySelector('title').textContent = 'bwa-mem2 on Biowulf';
bwa-mem2 on Biowulf


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



Bwa-mem2 is the next version of the bwa-mem algorithm in bwa. It produces alignment identical to bwa and is ~1.3-3.1x faster depending on the use-case, dataset and the running machine.


Reference

Vasimuddin Md, Sanchit Misra, Heng Li, Srinivas Aluru. Efficient Architecture-Aware Acceleration of BWA-MEM for Multicore Systems. IEEE Parallel and Distributed Processing Symposium (IPDPS), 2019. 
Documentation
* <https://github.com/bwa-mem2/bwa-mem2/tree/master>


Important Notes
* Module Name: bwa-mem2 (see [the modules page](/apps/modules.html) for more information)


* Multi-threaded using '-t' flag
* Index files in /fdb/bwa-mem2


Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive --mem=100g --cpus-per-task=4**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load bwa-mem2**

[user@cn3144 ~]$ **bwa-mem2 mem -t $SLURM\_CPUS\_PER\_TASK /fdb/bwa-mem2/hg38/genome.fa f1.fastq.gz f2.fastq.gz > /data/$USER/Ourfile**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. bwa-mem2.sh). For example:



```

#!/bin/bash
set -e
module load bwa-mem2
cd /data/$USER
bwa-mem2 mem -t $SLURM_CPUS_PER_TASK /fdb/bwa-mem2/hg38/genome.fa f1.fastq.gz f2.fastq.gz > /data/$USER/Ourfile

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --mem=50g --cpus-per-task=4 bwa-mem2.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. bwa-mem2.swarm). For example:



```

cd dir1; bwa-mem2 mem -t $SLURM_CPUS_PER_TASK /fdb/bwa-mem2/hg38/genome.fa f1.fastq.gz f2.fastq.gz > /data/$USER/Ourfile
cd dir2; bwa-mem2 mem -t $SLURM_CPUS_PER_TASK /fdb/bwa-mem2/hg38/genome.fa f1.fastq.gz f2.fastq.gz > /data/$USER/Ourfileq
...
cd dir10; bwa-mem2 mem -t $SLURM_CPUS_PER_TASK /fdb/bwa-mem2/hg38/genome.fa f1.fastq.gz f2.fastq.gz > /data/$USER/Ourfileq

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f bwa-mem2.swarm -g 50 -t 4 --module bwa-mem2
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module bwa-mem2 Loads the bwa-mem2 module for each subjob in the swarm 
 | |
 | |
 | |










