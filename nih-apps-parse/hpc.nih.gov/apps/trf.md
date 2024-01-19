

document.querySelector('title').textContent = 'TRF on Biowulf';
TRF on Biowulf


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



Tandem Repeats Finder is a program to locate and display tandem repeats in DNA sequences.



### References:


* G. Benson,
"Tandem repeats finder: a program to analyze DNA sequences"
Nucleic Acids Research (1999)
Vol. 27, No. 2, pp. 573-580. doi:[10.1093/nar/27.2.573](https://doi.org/10.1093/nar/27.2.573)


Documentation
* [TRF Main Site](https://tandem.bu.edu/trf/trf.html)
* [TRF Program Documentation](https://tandem.bu.edu/trf/trf.unix.help.html)


Important Notes
* Module Name: trf (see [the modules page](/apps/modules.html) for more information)



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem 5g**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load trf**
[+] Loading trf, version 4.09...
[user@cn3144 ~]$ **trf \
 /fdb/igenomes/Homo\_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa \
 2 `# match weight` \
 7 `# mismatch penalty` \
 7 `# indel penalty` \
 80 `# match probability` \
 10 `# indel probability` \
 50 `# minimum alignment score to report` \
 500 `# maximum period size to report` \
 -f `# record flanking sequence` \
 -d `# produce detailed data file` \
 -m `# generate a masked sequence file` \
 -l 6 `# longest tandem repeat array expected (in Mbp)`**

Tandem Repeats Finder, Version 4.09
Copyright (C) Dr. Gary Benson 1999-2012. All rights reserved.

Loading sequence...
Allocating Memory...
Initializing data structures...
Computing TR Model Statistics...
Scanning Sequence 1...
Freeing Memory...
Resolving output...
Done.
Loading sequence file...
Allocating Memory...
Initializing data structures...
Computing TR Model Statistics...
...

....................

.....................

.....................
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. trf.sh). For example:



```

#!/bin/sh

module load trf || exit 1

trf \
 /fdb/igenomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa \
 2 `# match weight` \
 7 `# mismatch penalty` \
 7 `# indel penalty` \
 80 `# match probability` \
 10 `# indel probability` \
 50 `# minimum alignment score to report` \
 500 `# maximum period size to report` \
 -f `# record flanking sequence` \
 -d `# produce detailed data file` \
 -m `# generate a masked sequence file` \
 -l 6 `# longest tandem repeat array expected (in Mbp)`

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --mem 5g --time 4:00:00 trf.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. trf.swarm). For example:



```

cd sample1 && trf sample1.fa 2 7 7 80 10 50 500 -f -d -m
cd sample2 && trf sample2.fa 2 7 7 80 10 50 500 -f -d -m
cd sample3 && trf sample3.fa 2 7 7 80 10 50 500 -f -d -m
cd sample4 && trf sample4.fa 2 7 7 80 10 50 500 -f -d -m

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f trf.swarm [-g #] [-t #] --module trf
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module trf Loads the trf module for each subjob in the swarm 
 | |
 | |
 | |








