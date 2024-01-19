

document.querySelector('title').textContent = 'Spades on HPC';
Spades on HPC


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

  SPAdes St. Petersburg genome assembler is intended for both standard isolate 
 and single-cell MDA bacteria assemblies. 
 ### References:


* If you use SPAdes in your research, please include [Nurk, 
 Bankevich et al., 2013](http://link.springer.com/chapter/10.1007/978-3-642-37195-0_13) in your reference list. You can also add [Bankevich, 
 Nurk et al., 2012](http://online.liebertpub.com/doi/abs/10.1089/cmb.2012.0021) instead.


Documentation
* <http://spades.bioinf.spbau.ru/release3.6.1/manual.html#sec4>



Important Notes
* Module Name: spades (see [the modules 
 page](/apps/modules.html) for more information)
* Multithreaded
* Example files in /usr/local/apps/spades/3.11.0/share/spades/





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

[user@cn3144 ~]$ **module load spades**[user@cn3144 ~]$ **spades.py -t $SLURM\_CPUS\_PER\_TASK -o outfile -1 infile\_1.fq.gz -2 infile\_2.fq.gz**
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```




Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. spades.sh). For example:



```

#!/bin/bash
set -e
module load spades
spades.py -t $SLURM_CPUS_PER_TASK -o outfile -1 infile_1.fq.gz -2 infile_2.fq.gz
```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=4 spades.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. spades.swarm). For example:



```

cd dir1; spades.py -t $SLURM_CPUS_PER_TASK -o outfile -1 infile_1.fq.gz -2 infile_2.fq.gz
cd dir2; spades.py -t $SLURM_CPUS_PER_TASK -o outfile -1 infile_1.fq.gz -2 infile_2.fq.gz
cd dir3; spades.py -t $SLURM_CPUS_PER_TASK -o outfile -1 infile_1.fq.gz -2 infile_2.fq.gz

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f spades.swarm -t 4 --module spades
```

where
 

|  |  |
| --- | --- |
| -t *#* | Number of threads/CPUs required for each process (1 line 
 in the swarm command file).  |
| --module  | Loads the module for each subjob in the swarm  |






