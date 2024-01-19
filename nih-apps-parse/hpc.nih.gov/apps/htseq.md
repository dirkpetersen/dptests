

document.querySelector('title').textContent = 'HTSeq on Biowulf';
HTSeq on Biowulf


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

  HTSeq is a Python package that provides infrastructure to process data from 
 high-throughput sequencing assays. It is developed by [Simon 
 Anders](http://www.embl.de/research/units/genome_biology/huber/members/index.php?s_personId=6001) at EMBL Heidelberg.



### References:

 * <https://academic.oup.com/bioinformatics/article/31/2/166/2366196>


Documentation * <http://www-huber.embl.de/users/anders/HTSeq/doc/tour.html#tour>


Important Notes
* Module Name: htseq (see [the modules 
 page](/apps/modules.html) for more information) 
 * Example files in <http://www-huber.embl.de/users/anders/HTSeq/HTSeq_example_data.tgz>



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

[user@cn3144 ~]$ **module load htseq**
[user@cn3144 ~]$ **htseq-count input.sam input.gff**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit). Create a batch input file (e.g. htseq.sh). For example:



```

#!/bin/bash
set -e
module load htseq
htseq-count input.sam input.gff
```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--mem=#] htseq.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. htseq.swarm). For example:



```

cd dir1; htseq-count input.sam input.gff
cd dir2; htseq-count input.sam input.gff
cd dir3; htseq-count input.sam input.gff
cd dir4; htseq-count input.sam input.gff

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f htseq.swarm [-t #] --module htseq
```

where
 

|  |  |  |  |
| --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in 
 the swarm command file) 
 | --module htseq Loads the htseq module for each subjob in the swarm 
  | |
 | |








