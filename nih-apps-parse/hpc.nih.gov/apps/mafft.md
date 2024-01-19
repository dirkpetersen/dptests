

document.querySelector('title').textContent = 'Mafft on HPC';
Mafft on HPC


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


MAFFT is a multiple sequence alignment program for unix-like operating 
 systems. It offers a range of multiple alignment methods, L-INS-i (accurate' 
 for alignment of < ~200 sequences), FFT-NS-2 (fast; for alignment of 
 < ~30000 sequences), etc


### References:


* [https://academic.oup.com/nar/article/30/14/3059/2904316/](https://mafft.cbrc.jp/alignment/server/jump.html?https://academic.oup.com/nar/article/30/14/3059/2904316/)


Documentation
* <https://mafft.cbrc.jp/alignment/software/>



Important Notes
* Module Name: mafft (see [the modules 
 page](/apps/modules.html) for more information)
* Multithreaded
* Example files in MAFFT\_TEST\_DATA





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

[user@cn3144 ~]$ **module load mafft**[user@cn3144 ~]$ **mafft --thread $SLURM\_CPUS\_PER\_TASK $MAFFT\_TEST\_DATA/ex1 > outfile**
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```




Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. batch.sh). For example:



```

#!/bin/bash
set -e
module load mafft
mafft --thread $SLURM_CPUS_PER_TASK ex1 > outfile

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=4 [--mem=#] batch.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. job.swarm). For example:



```

cd dir1; mafft --thread $SLURM_CPUS_PER_TASK ex1 > outfile
cd dir2; mafft --thread $SLURM_CPUS_PER_TASK ex1 > outfile
cd dir3; mafft --thread $SLURM_CPUS_PER_TASK ex1 > outfile
cd dir4; mafft --thread $SLURM_CPUS_PER_TASK ex1 > outfile

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f job.swarm [-g #] -t 4 --module mafft
```

where
 

|  |  |
| --- | --- |
| -g *#*  | Number of Gigabytes of memory required for each process (1 line in the swarm command file)
  |
| -t *#* | Number of threads/CPUs required for each process (1 line in the swarm command file).
  |
| --module  | Loads the module for each subjob in the swarm  |




