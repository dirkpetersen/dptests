

document.querySelector('title').textContent = 'raxml-ng on Biowulf';
raxml-ng on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
[Chaining long jobs](#chain)
[Swarm of jobs](#swarm) 
 |



RAxML-NG is a phylogenetic tree inference tool which uses maximum-likelihood (ML) optimality criterion. Its search heuristic is based on iteratively 
performing a series of Subtree Pruning and Regrafting (SPR) moves, which allows to quickly navigate to the best-known ML tree. RAxML-NG is a 
successor of RAxML (Stamatakis 2014) and leverages the highly optimized likelihood computation implemented in libpll (Flouri et al. 2014).

RAxML-NG offers improvements in speed, flexibility and user-friendliness over the previous RAxML versions. It also implements some of the features 
previously available in ExaML (Kozlov et al. 2015), including checkpointing and efficient load balancing for partitioned alignments (Kobert et al. 2014).



Reference

[RAxML-NG: a fast, scalable and user-friendly tool for maximum likelihood phylogenetic inference](https://academic.oup.com/bioinformatics/article/35/21/4453/5487384) 
Alexey M Kozlov, Diego Darriba, Tomás Flouri, Benoit Morel, Alexandros Stamatakis . Bioinformatics, Volume 35, Issue 21, 1 November 2019, Pages 4453-4455, https://doi.org/10.1093/bioinformatics/btz305


Documentation
* [raxml-ng on github](https://github.com/amkozlov/raxml-ng)
* [raxml github wiki](https://github.com/amkozlov/raxml-ng/wiki)


Important Notes
* Module Name: raxml-ng (see [the modules page](/apps/modules.html) for more information)
* Multithreaded, MPI
* Example files in /usr/local/apps/raxml-ng/examples



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

[user@cn3144 ~]$ **module load raxml-ng**

[user@cn3144 ~]$  **cp /usr/local/apps/raxml-ng/examples/myoglobin61.phy .**

[user@cn3144 ~]$  **raxml-ng --msa myoglobin597.phy --model GTR+G --threads $SLURM\_CPUS\_PER\_TASK**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. raxml-ng.sh).



**Threaded version of raxml-ng**

```

#!/bin/bash
set -e
module load raxml-ng

raxml-ng --msa myoglobin597.phy --model GTR+G --threads $SLURM_CPUS_PER_TASK

```

Submit this job using the SLURM [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] [--time=DD-HH:MM:SS] raxml-ng.sh
```


**MPI version of raxml-ng**

In our tests, the threaded version ran faster than the MPI version. However, for larger input datasets the MPI version might perform better. Please run your own benchmark comparisons for your data, and let us know of any interesting results.



Here is sample batch script for an MPI-only job (but on a single node):



```

#!/bin/bash

module load raxml-ng/1.0.0-mpi
mpirun -n $SLURM_NTASKS raxml-ng-mpi --msa data.phy --model GTR+G --threads 1

```


Submit with:

```

sbatch --ntasks=# --ntasks-per-core=1 [--time=DD-HH:MM:SS] jobscript

```

Make sure to consult our documentation on [efficient use of the multinode partition](https://hpc.nih.gov/policies/multinode.html) if you want to run raxml-ng-mpi over multiple nodes. Before you use the multinode parition make sure your dataset could benefit from such a large number of resources.


In most cases where you want to use MPI, the hybrid MPI/pthreads setup is more efficient. Please consult [this page](https://github.com/amkozlov/raxml-ng/wiki/Parallelization#mpi-and-hybrid-mpipthreads) for how to set these up and reach out to us if needed.



Chaining long raxml-ng jobs

Checkpointing is built in to raxml-ng. Thus, if you have a longrunning raxml-ng job that is terminated because you did not specify a long enough walltime, or if your 
job is going to require more than the Biowulf 10-day max walltime, you can resubmit the job and it will restart from the last checkpoint. 
See [the advanced tutorial for details](https://github.com/amkozlov/raxml-ng/wiki/Advanced-Tutorial).

You can take advantage of this feature to run a chain of jobs, each of which will pick up where the previous one terminated, using [job dependencies](https://hpc.nih.gov/docs/userguide.html#depend). For example:


```

biowulf% **sbatch --cpus-per-task=16 --mem=20g --time=8-00:00:00 myjobscript**
1111
biowulf% **sbatch --depend=afterany:1111 --cpus-per-task=16 --mem=20g --time=8-00:00:00 myjobscript**
2222
biowulf% **sbatch --depend=afterany:2222 --cpus-per-task=16 --mem=20g --time=8-00:00:00 myjobscript**
3333

```

Each of these jobs will run for 8 days. When job 1111 terminates, job 2222 will start up from the last checkpoint file, and likewise for job 3333. The 3 jobs will utilize a total walltime of 24 days. 



Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. raxml-ng.swarm). For example:



```

raxml-ng --msa file1.phy --model GTR+G --threads $SLURM_CPUS_PER_TASK
raxml-ng --msa file2.phy --model GTR+G --threads $SLURM_CPUS_PER_TASK
raxml-ng --msa file3.phy --model GTR+G --threads $SLURM_CPUS_PER_TASK
[...]

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f raxml-ng.swarm -g 20  -t 8  --module raxml-ng
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g 20 20 Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t 8 8 threads/CPUs required for each process (1 line in the swarm command file).
 | --module raxml-ng Loads the raxml-ng module for each subjob in the swarm 
 | |
 | |
 | |


























