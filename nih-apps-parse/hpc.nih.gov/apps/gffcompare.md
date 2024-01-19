

document.querySelector('title').textContent = 'gffcompare on Biowulf';
gffcompare on Biowulf


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



The program gffcompare can be used to compare, merge, annotate and estimate accuracy of one or more GFF files (the “query” files), when compared with a reference annotation (also provided as GFF).

This program is based on the CuffCompare utility which is part of the Cufflinks/Tuxedo suite, so the various usage options and output files as documented in the CuffCompare manual apply to the gffcompare program as well.




Documentation
* [gffcompare Main Site](http://ccb.jhu.edu/software/stringtie/gffcompare.shtml)
* [gffcompare on GitHub](https://github.com/gpertea/gffcompare)


Important Notes
* Module Name: gffcompare (see [the modules page](/apps/modules.html) for more information)



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

[user@cn3144 ~]$ **module load gffcompare**

[user@cn3144 ~]$ **gffcompare -R -r mm10.gff -o cuffcmp cufflinks\_asm.gtf**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. gffcompare.sh). For example:



```

#!/bin/bash
set -e
module load gffcompare
gffcompare -R -r mm10.gff -o cuffcmp cufflinks_asm.gtf

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] gffcompare.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. gffcompare.swarm). For example:



```

gffcompare -R -r mm10.gff -o cuffcmp1 cufflinks_asm1.gtf
gffcompare -R -r mm10.gff -o cuffcmp2 cufflinks_asm2.gtf
gffcompare -R -r mm10.gff -o cuffcmp3 cufflinks_asm3.gtf
gffcompare -R -r mm10.gff -o cuffcmp4 cufflinks_asm4.gtf

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f gffcomapre.swarm [-g #] [-t #] --module gffcompare
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module gffcompare Loads the gffcompare module for each subjob in the swarm 
 | |
 | |
 | |








