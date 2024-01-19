

document.querySelector('title').textContent = 'Pear on HPC';
Pear on HPC


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


 **PEAR** is an ultrafast, memory-efficient and highly accurate 
 pair-end read merger. It is fully parallelized and can run with as low as 
 just a few kilobytes of memory. 


PEAR evaluates all possible paired-end read overlaps and without requiring 
 the target fragment size as input. In addition, it implements a statistical 
 test for minimizing false-positive results. 


### 


Documentation
* <https://sco.h-its.org/exelixis/web/software/pear/doc.html>



Important Notes
* Module Name: pear (see [the modules 
 page](/apps/modules.html) for more information)
* Multithreaded





Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --cpus-per-task=4**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load pear**
[user@cn3144 ~]$ **pear -f input\_R1.fastq.gz -r input\_R2.fastq.gz -o output -j $SLURM\_CPUS\_PER\_TASK**
 ____  _____    _    ____
|  _ \| ____|  / \  |  _ \
| |_) |  _|   / _ \ | |_) |
|  __/| |___ / ___ \|  _ <
|_|   |_____/_/   \_\_| \_\

PEAR v0.9.11 [Nov 5, 2017]

Citation - PEAR: a fast and accurate Illumina Paired-End reAd mergeR
Zhang et al (2014) Bioinformatics 30(5): 614-620 | doi:10.1093/bioinformatics/btt593

Forward reads file.................: input_R1.fastq.gz
Reverse reads file.................: input_R2.fastq.gz
PHRED..............................: 33
Using empirical frequencies........: YES
Statistical method.................: OES
Maximum assembly length............: 999999
Minimum assembly length............: 50
p-value............................: 0.010000
Quality score threshold (trimming).: 0
Minimum read size after trimming...: 1
Maximal ratio of uncalled bases....: 1.000000
Minimum overlap....................: 10
Scoring method.....................: Scaled score
Threads............................: 4

Allocating memory..................: 200,000,000 bytes
Computing empirical frequencies....: DONE
  A: 0.256699
  C: 0.237672
  G: 0.250363
  T: 0.255267
  693857 uncalled bases
Assemblying reads: 100%

Assembled reads ...................: 9,010,229 / 10,000,000 (90.102%)
Discarded reads ...................: 0 / 10,000,000 (0.000%)
Not assembled reads ...............: 989,771 / 10,000,000 (9.898%)
Assembled reads file...............: output.assembled.fastq
Discarded reads file...............: output.discarded.fastq
Unassembled forward reads file.....: output.unassembled.forward.fastq
Unassembled reverse reads file.....: output.fastq.unassembled.reverse.fastq

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```




Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. batch.sh). For example:



```

#!/bin/bash
set -e
module load pear
pear -f input1.fastq.gz -r input2.fastq.gz -o output -j $SLURM_CPUS_PER_TASK
```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=4 batch.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. job.swarm). For example:



```

cd dir1; pear -f S1_R1.fastq.gz -r S1_R2.fastq.gz -o S1 -j $SLURM_CPUS_PER_TASK
cd dir2; pear -f S2_R1.fastq.gz -r S2_R2.fastq.gz -o S2 -j $SLURM_CPUS_PER_TASK
cd dir3; pear -f S3_R1.fastq.gz -r S3_R2.fastq.gz -o S3 -j $SLURM_CPUS_PER_TASK

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f job.swarm -t 4 --module pear
```

where
 

|  |  |
| --- | --- |
| -t *#* | Number of threads/CPUs required for each process (1 line 
 in the swarm command file).  |
| --module  | Loads the module for each subjob in the swarm  |




