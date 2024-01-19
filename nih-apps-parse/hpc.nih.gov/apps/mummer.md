

document.querySelector('title').textContent = 'MUMmer on Biowulf';
MUMmer on Biowulf


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



MUMmer is a system for rapidly aligning entire genomes. 



**Reference:**  

[Versatile and open software for comparing large genomes](https://www.ncbi.nlm.nih.gov/pubmed/14759262)
S. Kurtz, A. Phillippy, A.L. Delcher, M. Smoot, M. Shumway, C. Antonescu, and S.L. Salzberg.
Genome Biology (2004), 5:R12.

Documentation
* [MUMmer Manual](https://github.com/mummer4/mummer/blob/master/README.md)


Important Notes
* Module Name: mummer (see [the modules page](/apps/modules.html) for more information)
* Multithreaded.
* mummerplot fails with any gnuplot later than 4.6.5. So you should 'module load gnuplot/4.6.5' before running mummerplot.



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

[user@cn3144 ~]$ **module load mummer**

[user@cn3144 ~]$ **mummer /fdb/genome/hg19/chr\_all.fa /fdb/genome/hg19/chrX.fa**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. MUMmer.sh). For example:



```

#!/bin/bash
set -e
module load mummer
mummer /fdb/genome/hg19/chr_all.fa /fdb/genome/hg19/chrX.fa

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --mem=10g MUMmer.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. MUMmer.swarm). For example:



```

mummer    /fdb/genome/hg19/chr_all.fa     file1.fa
mummer    /fdb/genome/hg19/chr_all.fa     file2.fa
mummer    /fdb/genome/hg19/chr_all.fa     file3.fa

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f MUMmer.swarm -g 1-0 [-t #] --module MUMmer
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module MUMmer Loads the MUMmer module for each subjob in the swarm 
 | |
 | |
 | |








