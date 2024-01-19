

document.querySelector('title').textContent = 'bcl2fastq on Biowulf';
bcl2fastq on Biowulf


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


The primary output of Illumina sequencing instruments are per-cycle base call 
files in BCL format. Bcl2fastq converts BCL files to fastq files used by most
downstream software. It also separates reads into individual fastq files based
on their barcode (demultiplexing), does adapter masking/trimming and moves unique
molecular identifier (UMI) bases from the read to the fastq header. Bcl2fastq
operates on Illumina run folders.


bcl2fastq creates a number of files summarizing statistics of the conversion
in the InterOp folder of the run folder. They can be examined with the Illumina
Sequence Analysis Viewer.


Documentation
* [User guide for v2.20 [PDF]](https://support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/bcl2fastq/bcl2fastq2-v2-20-software-guide-15051736-03.pdf)* [User guide for version 2.19 [PDF]](https://support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/bcl2fastq/bcl2fastq2_guide_15051736_v2.pdf)


Important Notes
* Module Name: bcl2fastq (see [the modules page](/apps/modules.html) for more information)


bcl2fastq is a multithreaded application. By default it will run as many
threads as there are CPUs for the conversion/demultiplexing plus additional
threads for reading/writing data. **You must therefore limit the threads to the number of threads
allocated to your job or allocate nodes exclusively.**


 The relevant options for limiting the number of threads are




```

-r [ --loading-threads ] arg (=4)     number of threads used for loading BCL data
-p [ --processing-threads ] arg       number of threads used for processing demultiplexed data
-w [ --writing-threads ] arg (=4)     number of threads used for writing FASTQ data

```

Version 2.17 has an additional setting for the numer demultiplexing threads which is
not present in later versions:



```

-d [ --demultiplexing-threads ] arg   number of threads used for demultiplexing

```

For a job allocated 16 CPUs these options should be set to



```

bcl2fastq -r 4 -w 4 -p 14 ...

```

This leads to a nominal overload but according to testing by Illumina
this should yield optimal throughput.



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --cpus-per-task=16**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load bcl2fastq**
[user@cn3144 ~]$ **bcl2fastq --help**
BCL to FASTQ file converter
bcl2fastq v2.20.0.422
Copyright (c) 2007-2017 Illumina, Inc.

Usage:
      bcl2fastq [options]
[...snip...]

[user@cn3144 ~]$ **bcl2fastq -r 4 -w 4 -p 14 \
 --runfolder-dir /path/to/your/run/folder/160601\_mach1\_0023 ...**

[user@cn3144 ~]$

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. bcl2fastq.sh), which uses the input file 'bcl2fastq.in'. For example:



```

#! /bin/bash

module load bcl2fastq/2.20.0 || exit 1
bcl2fastq --runfolder-dir /path/to/your/run/folder/160601_mach1_0023 \
          --output-dir ./160601_mach1_0023 \
          -r 4 -w 4 -p 14 \
          --barcode-mismatches 0

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=16 bcl2fastq.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. bcl2fastq.swarm). For example:



```

bcl2fastq --runfolder-dir /path/to/your/run/folder/160601_mach1_0023 \
          --output-dir ./160601_mach1_0023 \
          -r 4 -w 4 -p 14
bcl2fastq --runfolder-dir /path/to/your/run/folder/160602_mach1_0024 \
          --output-dir ./160602_mach1_0024 \
          -r 4 -w 4 -p 14
bcl2fastq --runfolder-dir /path/to/your/run/folder/160603_mach1_0025 \
          --output-dir ./160603_mach1_0025 \
          -r 4 -w 4 -p 14

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f bcl2fastq.swarm -t 16 --module bcl2fastq
```

where


|  |  |  |  |
| --- | --- | --- | --- |
| -t #  Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module bcl2fastq  Loads the bcl2fastq module for each subjob in the swarm 
 | |
 | |








