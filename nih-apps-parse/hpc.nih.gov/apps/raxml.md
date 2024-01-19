

document.querySelector('title').textContent = 'raxml on Biowulf';
raxml on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
[Swarm of jobs](#swarm) 
[MPI Parallel job](#MPI)
 |



RAxML (Randomized Axelerated Maximum Likelihood) is a program for sequential and parallel
Maximum Likelihood based inference of large phylogenetic trees. It can also be used for post-
analyses of sets of phylogenetic trees, analyses of alignments and, evolutionary placement of short
reads.
It has originally been derived from fastDNAml which in turn was derived from Joe Felsentein's
dnaml which is part of the PHYLIP package. 
[[RAxML website](http://sco.h-its.org/exelixis/web/software/raxml/index.html)]


There are several different RAxML executables:


|  |  |  |  |  |  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| | raxmlHPC  sequential version. Intended for small to medium datasets and for initial experiments to determine appropriate search parameters. This is not built with SSE3, and for performance reasons one of the other versions is preferred. 
| raxmlHPC-SSE3  sequential version built with SSE3. Intended for small to medium datasets and for initial experiments to determine appropriate search parameters.
| raxmlHPC-PTHREADS  Can run multiple threads on multiple cores of a single node. Works well for longer alignments..
| raxmlHPC-PTHREADS-SSE3  Can run multiple threads on multiple cores of a single node, built with SSE3. Works well for longer alignments..
| raxmlHPC-MPI  Can run multiple MPI processes on multiple cores of multiple nodes. Intended for executing really large production runs 
(i.e. 100 or 1,000 bootstraps). It has been designed to do multiple
inferences or rapid/standard BS (bootstrap) searches in parallel! 
For all remaining options, the usage of this type of coarse-grained parallelism does not make much
sense!
 | |
 | |
 | |
 | |
 | |
 |





Documentation
* [raxml Manual](http://sco.h-its.org/exelixis/php/countManualNew.php)
* [RAxML memory calculator](https://cme.h-its.org/exelixis/web/software/raxml/).


Important Notes
* Module Name: raxml (see [the modules page](/apps/modules.html) for more information)
* read the section in the [user manual](http://sco.h-its.org/exelixis/php/countManualNew.php) on how many threads/cores to use. 
Rough rule of thumb: 1 thread/core per 500 DNA site patterns.
* environment variables set 
	+ raxml\_HOME* The test datasets for RAxML are available in /usr/local/apps/raxml/test-data.



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

[user@cn3144 ~]$ **module load raxml**

[user@cn3144 ~]$ **raxmlHPC -m BINGAMMA -p 12345 -s binary.phy -n T1** 
IMPORTANT WARNING: Sequences t2 and t3 are exactly identical
IMPORTANT WARNING: Sequences t2 and t4 are exactly identical
IMPORTANT WARNING
Found 2 sequences that are exactly identical to other sequences in the alignment.
Normally they should be excluded from the analysis.

Just in case you might need it, an alignment file with 
sequence duplicates removed is printed to file binary.phy.reduced

[...]
Starting final GAMMA-based thorough Optimization on tree 0 likelihood -119.520001 .... 
Final GAMMA-based Score of best tree -119.520001

Program execution info written to /spin1/users/user/raxml/RAxML_info.T1
Best-scoring ML tree written to: /spin1/users/user/raxml/RAxML_bestTree.T1
Overall execution time: 0.496108 secs or 0.000138 hours or 0.000006 days

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. raxml.sh). For example:



```

#!/bin/bash
set -e
module load raxml
raxmlHPC-PTHREADS -m BINGAMMA -p 12345 -s binary.phy -n T3  -T $SLURM_CPUS_PER_TASK

```

The -T $SLURM\_CPUS\_PER\_TASK flag specifies the number of threads to run. This will automatically be the same as the
number of CPUs you allocate using the sbatch command with --cpus-per-task, as in the example above. As per the Raxml
manual, it is very important to run the same number of threads as CPUs. 

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=# [--mem=#] raxml.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. raxml.swarm). For example:



```

raxmlHPC -m BINGAMMA -p 12345 -s file1.phy -n T3  
raxmlHPC -m BINGAMMA -p 12345 -s file2.phy -n T3  
raxmlHPC -m BINGAMMA -p 12345 -s file3.phy -n T3  
[...]     

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f raxml.swarm [-g #] --module raxml
```

where


|  |  |  |  |
| --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | --module raxml Loads the raxml module for each subjob in the swarm 
 | |
 | |


If you are using threads (the 'raxmlHPC-PTHREADS' executable and -T $SLURM\_CPUS\_PER\_TASK), be sure to specify the number of threads to swarm with:

```

$ swarm -f cmdfile -t # --module raxml

```


 MPI RAxML batch job

An MPI job can run in parallel across multiple nodes. It can also run on multiple CPUs of a single node, similar to the threaded version.
Set up a batch script along the following lines:

```

#!/bin/bash
#PBS -N raxml

cd /data/user/raxml

module load raxml

echo "Running $SLURM_NTASKS MPI processes "

mpirun -np $SLURM_NTASKS raxmlHPC-MPI -f d -m GTRCAt -p 12345 -s ex_al -N 100 -n MultipleOriginal

```

Submit the job with:

```

sbatch --ntasks=8 --ntasks-per-core=1 jobscript

```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| --ntasks=8  tells slurm you need to run 8 MPI processes (up to 28 tasks can be submitted to a single node in the norm partition, and for any more than 28 you should submit to the multinode partition)
| --ntasks-per-core=1  run only one task on a physical core (usually recommended for MPI jobs)
| --np $SLURM\_NTASKS  tells mpirun to run $SLURM\_NTASKS processes, which is set to 8 via the sbatch command line. 
 | |
 | |
 | |


You will need to experiment with the number of tasks to find the optimal settings. The batch script can remain the same, and you can modify the sbatch command line to try different values. The slurm-####.out file will report the number of tasks, and you can use [jobhist](/docs/biowulf_tools.html#jobhist) to see how long a job took to complete. Please send us the results of your experiment -- i.e. the time taken for various values of ntasks -- to help in making recommendations to other users. 













