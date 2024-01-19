

document.querySelector('title').textContent = 'alleleCount on Biowulf';
alleleCount on Biowulf


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



AlleleCount provides support code for NGS copy number algorithms. It exists primarily to prevent code duplication between some other projects, specifically AscatNGS and Battenburg.
The project contains 2 equivalent implementations of allele counting code in perl and C for BAM/CRAM processing.



Documentation
* [alleleCount website](https://github.com/cancerit/alleleCount)


Important Notes
* Module Name: alleleCount (see [the modules page](/apps/modules.html) for more information)
* Singletheaded
* Environment variable set, $ALLELECOUNT\_TEST\_DATA



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

[user@cn3144 ~]$ **module load alleleCount**

[user@cn3144 ~]$ **cp $ALLELECOUNT\_TEST\_DATA/loci\_22.txt .**

[user@cn3144 ~]$ **cp $ALLELECOUNT\_TEST\_DATA/test.bam\* .**

[user@cn3144 ~]$ **alleleCounter -l loci\_22.txt -b test.bam -o test.out**
Reading locis
Done reading locis

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. alleleCount.sh). For example:



```

#!/bin/bash
set -e
cd /data/$USER/mydir
module load alleleCount
alleleCounter -l loci.txt -b input.bam -o output.txt


```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--mem=#] alleleCount.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. alleleCount.swarm). For example:



```

alleleCounter -l loci1.txt -b input.bam -o output1.txt
alleleCounter -l loci2.txt -b input.bam -o output2.txt
[...]

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f alleleCount.swarm [-g #] --module alleleCount
```

where


|  |  |  |  |
| --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | --module alleleCount Loads the alleleCount module for each subjob in the swarm 
 | |
 | |








