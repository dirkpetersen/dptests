

document.querySelector('title').textContent = 'porechop on Biowulf';
porechop on Biowulf


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


Porechop finds and removes adapters from Oxford Nanopore reads. Adapters on
the ends of reads are trimmed off. Reads with adapters in its middle are
treated as chimeric and split into separate reads. Porechop performs thorough
alignments to effectively find adapters, even at low sequence identity.


Reads barcoded with some barcoding kits can also be demultiplexed.


Documentation
* porechop [on GitHub/rrwick/Porechop](https://github.com/rrwick/Porechop)


Important Notes
* Module Name: porechop (see [the modules page](/apps/modules.html) for more information)
* porechop can use multiple threads to increase processing speed



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

[user@cn3144 ~]$ **module load porechop**
[user@cn3144 ~]$ **porechop -i input\_reads.fastq.gz -o output\_reads.fastq.gz**
...
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. porechop.sh), which uses the input file 'porechop.in'. For example:



```

#!/bin/bash
module load porechop/0.2.3 || exit 1
porechop -i input_reads.fastq.gz -b output_dir --threads=$SLURM_CPUS_PER_TASK

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=2 [--mem=#] porechop.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. porechop.swarm). For example:



```

porechop -i input_reads1.fastq.gz -o output_reads1.fastq.gz --verbosity 2 --threads $SLURM_CPUS_PER_TASK
porechop -i input_reads2.fastq.gz -o output_reads2.fastq.gz --verbosity 2 --threads $SLURM_CPUS_PER_TASK
porechop -i input_reads3.fastq.gz -o output_reads3.fastq.gz --verbosity 2 --threads $SLURM_CPUS_PER_TASK

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f porechop.swarm [-g #] -t 6 --module porechop
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g #  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t #  Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module porechop  Loads the porechop module for each subjob in the swarm 
 | |
 | |
 | |








