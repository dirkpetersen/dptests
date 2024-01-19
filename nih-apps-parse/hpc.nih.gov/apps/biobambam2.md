

document.querySelector('title').textContent = 'BIOBAMBAM2 on Biowulf';
BIOBAMBAM2 on Biowulf


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



Biobambam2 is a toolkit that contains tools for early stage alignment file processing. The following tools are provided: 
* bamsormadup: parallel sorting and duplicate marking
* bamcollate2: reads BAM and writes BAM reordered such that alignment or collated by query name
* bammarkduplicates: reads BAM and writes BAM with duplicate alignments marked using the BAM flags field
* bammaskflags: reads BAM and writes BAM while masking (removing) bits from the flags column
* bamrecompress: reads BAM and writes BAM with a defined compression setting. This tool is capable of multi-threading.
* bamsort: reads BAM and writes BAM resorted by coordinates or query name
* bamtofastq: reads BAM and writes FastQ; output can be collated or uncollated by query name





Documentation
* [biobambam2 Main Site](https://gitlab.com/german.tischler/biobambam2)


Important Notes
* Module Name: biobambam2 (see [the modules page](/apps/modules.html) for more information)



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session with bamsort example (user input in **bold**):



```

[user@biowulf]$ **sinteractive**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load biobambam2**

[user@cn3144 ~]$ **bamsort SO=queryname <wgEncodeUwRepliSeqBg02esG1bAlnRep1.bam >out.bam**
[V] Reading alignments from source.
[V] 1M
[V] read 1784867 alignments
[V] producing sorted output
[V]1
[V] wrote 1784867 alignments

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. biobambam2.sh). For example:



```

#!/bin/bash
set -e
module load biobambam2
cd /data/username
bamsort SO=queryname <wgEncodeUwRepliSeqBg02esG1bAlnRep1.bam >out.bam

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] biobambam2.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. biobambam2.swarm). For example:



```

bamsort SO=queryname <input1.bam >output1.bam
bamsort SO=queryname <input2.bam >output2.bam
bamsort SO=queryname <input3.bam >output3.bam
bamsort SO=queryname <input4.bam >output4.bam

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f biobambam2.swarm [-g #] [-t #] --module biobambam2
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module biobambam2 Loads the biobambam2 module for each subjob in the swarm 
 | |
 | |
 | |








