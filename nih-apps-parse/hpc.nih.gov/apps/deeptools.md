

document.querySelector('title').textContent = 'Deeptools on Biowulf';
Deeptools on Biowulf


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

  deepTools is a suite of user-friendly tools for the visualization, quality 
 control and normalization of data from high-throughput DNA sequencing experiments. 
 


deepTools offers multiple methods for highly-customizable data visualization 
 that immensely aid hypothesis generation and data interpretation. It also 
 offers all the tools needed to create coverage files in standard bedGraph 
 and bigWig file formats allowing various normalization procedures and comparisons 
 between two files (for example, treatment and control).


 


### References:


* <https://academic.oup.com/nar/article/44/W1/W160/2499308>


Documentation * <http://deeptools.readthedocs.io/en/latest/>



Important Notes
* Module Name: deeptools (see [the 
 modules page](/apps/modules.html) for more information)
* Programs auto thread to all available processors. So make sure -p flag 
 is specified.





Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --cpus-per-task=4**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load deeptools**
[user@cn3144 ~]$ **bamCoverage -b file.bam -o outfile -of bigwig -p $SLURM\_CPUS\_PER\_TASK**
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```

$SLURM\_CPUS\_PER\_TASK will be automatically replaced by the number used in the sinteractive command (4 in this example)



Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. deeptools.sh). For example:



```

#!/bin/bash
set -e
module load deeptools
bamCoverage -b file.bam -o outfile -of bigwig -p $SLURM_CPUS_PER_TASK

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=4 --mem=10g deeptools.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. deeptools.swarm). For example:



```
cd dir1; bamCoverage -b file.bam -o outfile -of bigwig -p $SLURM_CPUS_PER_TASK
cd dir2; bamCoverage -b file.bam -o outfile -of bigwig -p $SLURM_CPUS_PER_TASK
cd dir3; bamCoverage -b file.bam -o outfile -of bigwig -p $SLURM_CPUS_PER_TASK
cd dir4; bamCoverage -b file.bam -o outfile -of bigwig -p $SLURM_CPUS_PER_TASK

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f deeptools.swarm g 10 t 4 --module deeptools
```

where
 

|  |  |
| --- | --- |
| -g *#*  | Number of Gigabytes of memory required for each process (1 line in the swarm command file)
  |
| -t *#* | Number of threads/CPUs required for each process (1 line in the swarm command file).
  |
| --module  | Loads the module for each subjob in the swarm  |




