

document.querySelector('title').textContent = 'FastTree on Biowulf';
FastTree on Biowulf


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



FastTree infers approximately-maximum-likelihood phylogenetic trees from alignments of nucleotide or protein sequences. FastTree can handle alignments with up to a million of sequences in a reasonable amount of time and memory. For large alignments, FastTree is 100-1,000 times faster than PhyML 3.0 or RAxML 7. 



### References:


* FastTree is open-source software developed by Price et al. [[FastTree paper](http://mbe.oxfordjournals.org/content/26/7/1641.full)]


Documentation
* [FastTree documentation](http://www.microbesonline.org/fasttree/)


Important Notes
* Module Name: FastTree (see [the modules page](/apps/modules.html) for more information)
* Singlethreaded (FastTree) and Multithreaded (FastTreeMP) versions are available.



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

[user@cn3144 ~]$ **module load FastTree**

[user@cn3144 ~]$ **FastTreeMP alignmentfile > treefile**
[....]

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. FastTree.sh). For example:



```

#!/bin/bash
set -e

module load FastTree

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

FastTreeMP alignmentfile > treefile

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=8 --mem=4g FastTree.sh
```


The command above allocates 8 CPUs for the job. Within the batch script, the variable $SLURM\_CPUS\_PER\_TASK is used to set the number of threads that FastTreeMP spawns. The number 
of threads should always be the same as the number of allocated CPUs, and you make sure of this by using $SLURM\_CPUS\_PER\_TASK within the script. 

If you need more than the default 4 GB of memory, add  --mem=8g to the sbatch command line.

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. FastTree.swarm). For example:



```

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK; FastTreeMP  alignment1 > tree1
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK; FastTreeMP  alignment2 > tree2
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK; FastTreeMP  alignment3 > tree3
[....]

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f FastTree.swarm [-g #] [-t #] --module FastTree
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module FastTree Loads the FastTree module for each subjob in the swarm 
 | |
 | |
 | |












