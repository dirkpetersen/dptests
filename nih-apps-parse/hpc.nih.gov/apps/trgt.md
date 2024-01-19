

document.querySelector('title').textContent = "trgt";
trgt on Biowulf


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



TRGT is a tool for targeted genotyping of tandem repeats from PacBio HiFi data. In addition to the basic size genotyping, TRGT profiles sequence composition, mosaicism, and CpG methylation of each analyzed repeat. TRGT comes with a companion tool TRVZ for visualization of reads overlapping the repeats.




Documentation
* [trgt Main Site](https://github.com/PacificBiosciences/trgt)


Important Notes
* Module Name: trgt (see [the modules page](/apps/modules.html) for more information)
 * Multithreaded. Set number of threads with trgt's --threads flag.
* trvz produces SVG images. You can use [hpcdrive](/docs/hpcdrive.html) to view them on your local workstation.
* Environment variables set 
	+ TRGT\_HOME* Example files in $TRGT\_HOME/example* Reference data in $TRGT\_HOME/repeats



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load trgt**

[user@cn3144 ~]$ **trgt --genome $TRGT\_HOME/example/reference.fasta --repeats $TRGT\_HOME/example/repeat.bed --reads $TRGT\_HOME/example/sample.bam --output-prefix sample**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. trgt.sh). For example:



```

#!/bin/bash
set -e
module load trgt
trgt --genome $TRGT_HOME/example/reference.fasta --repeats $TRGT_HOME/example/repeat.bed --reads $TRGT_HOME/example/sample.bam --output-prefix sample

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] trgt.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. trgt.swarm). For example:



```

trgt --genome $TRGT_HOME/example/reference.fasta --repeats $TRGT_HOME/example/repeat.bed --reads sample1.bam --output-prefix sample1
trgt --genome $TRGT_HOME/example/reference.fasta --repeats $TRGT_HOME/example/repeat.bed --reads sample2.bam --output-prefix sample2
trgt --genome $TRGT_HOME/example/reference.fasta --repeats $TRGT_HOME/example/repeat.bed --reads sample3.bam --output-prefix sample3
trgt --genome $TRGT_HOME/example/reference.fasta --repeats $TRGT_HOME/example/repeat.bed --reads sample4.bam --output-prefix sample4

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f trgt.swarm [-g #] [-t #] --module trgt
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module trgt Loads the trgt module for each subjob in the swarm 
 | |
 | |
 | |








