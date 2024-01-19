

document.querySelector('title').textContent = 'LASTZ on Biowulf';
LASTZ on Biowulf


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



LASTZ is a tool for (1) aligning two DNA sequences, and (2) inferring appropriate scoring parameters automatically. It is a drop-in replacement for BLASTZ, and is backward compatible with BLASTZ's command-line syntax.


LASTZ is designed to preprocess one sequence or set of sequences (which we collectively call the target) and then align several query sequences to it. The general flow of the program is like a pipeline: the output of one stage is the input to the next. The user can choose to skip most stages via command-line options; any stages that are skipped pass their input along to the next stage unchanged. Two of the stages, scoring inference and interpolation, are special in that they perform a miniature version of the pipeline within them.



### References:


* [Harris, R.S. (2007) Improved pairwise alignment of genomic DNA. Ph.D. Thesis, The Pennsylvania State University.](http://www.bx.psu.edu/~rsharris/rsharris_phd_thesis_2007.pdf)


Documentation
* [LASTZ Main Site](http://www.bx.psu.edu/~rsharris/lastz/)


Important Notes
* Module Name: LASTZ (see [the modules page](/apps/modules.html) for more information)
* Single-threaded
* Reference data in /data/genome/fasta/



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

[user@cn3144 ~]$ lastz hg18.chr4.fa galGal3.chr4.fa --notransition --step=20 --nogapped --format=maf > hg18_4_vs_galGal3_4.maf 

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. LASTZ.sh). For example:



```

#!/bin/bash
set -e
module load LASTZ
lastz hg18.chr4.fa galGal3.chr4.fa --notransition --step=20 --nogapped --format=maf > hg18_4_vs_galGal3_4.maf

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] LASTZ.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. LASTZ.swarm). For example:



```

lastz hg18.chr1.fa galGal3.chr1.fa --notransition --step=20 --nogapped --format=maf > hg18_1_vs_galGal3_1.maf
lastz hg18.chr2.fa galGal3.chr2.fa --notransition --step=20 --nogapped --format=maf > hg18_2_vs_galGal3_2.maf
lastz hg18.chr3.fa galGal3.chr3.fa --notransition --step=20 --nogapped --format=maf > hg18_3_vs_galGal3_3.maf
lastz hg18.chr4.fa galGal3.chr4.fa --notransition --step=20 --nogapped --format=maf > hg18_4_vs_galGal3_4.maf

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f LASTZ.swarm [-g #] [-t #] --module LASTZ
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module LASTZ Loads the LASTZ module for each subjob in the swarm 
 | |
 | |
 | |








