

document.querySelector('title').textContent = 'iqtree on Biowulf';
iqtree on Biowulf


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



Efficient software for phylogenomic inference

> 
> The IQ-TREE software was created as the successor of IQPNNI and TREE-PUZZLE (thus the name IQ-TREE).
> IQ-TREE was motivated by the rapid accumulation of phylogenomic data, leading to a need for efficient
> phylogenomic software that can handle a large amount of data and provide more complex models of
> sequence evolution. To this end, IQ-TREE can utilize multicore computers and distributed parallel
> computing to speed up the analysis. IQ-TREE automatically performs checkpointing to resume an
> interrupted analysis.
> 





### References:


* B.Q. Minh, H.A. Schmidt, O. Chernomor, D. Schrempf, M.D. Woodhams, A. von Haeseler, R. Lanfear (2020)
 [**IQ-TREE 2: New models and efficient methods for
 phylogenetic inference in the genomic era.**](https://doi.org/10.1093/molbev/msaa015)
*Mol. Biol. Evol., 37:1530-1534.*
* L.-T. Nguyen, H.A. Schmidt, A. von Haeseler, B.Q. Minh (2015)
 [**IQ-TREE: A fast and effective stochastic algorithm
 for estimating maximum likelihood phylogenies**](https://doi.org/10.1093/molbev/msu300)
*Mol. Biol. Evol., 32:268-274.*


Documentation
* [iqtree Main Site](http://www.iqtree.org/)
* [iqtree Github](https://github.com/iqtree/iqtree2)


Important Notes
* Module Name: iqtree (see [the modules page](/apps/modules.html) for more information)
* Multithreaded application
* Example files in $IQTREE\_TEST\_DATA



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=6g -c 8**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load iqtree**

[user@cn3144 ~]$ **cp $IQTREE\_TEST\_DATA/\* .**

[user@cn3144 ~]$ **iqtree2 -s example.phy -nt $SLURM\_CPUS\_PER\_TASK**
IQ-TREE multicore version 2.1.2 COVID-edition for Linux 64-bit built Oct 22 2020
Developed by Bui Quang Minh, James Barbetti, Nguyen Lam Tung,
Olga Chernomor, Heiko Schmidt, Dominik Schrempf, Michael Woodhams.

Host:    cn0852 (AVX512, FMA3, 377 GB RAM)
Command: iqtree2 -s example.phy -nt 8
Seed:    830607 (Using SPRNG - Scalable Parallel Random Number Generator)
Time:    Mon Oct 26 10:33:27 2020
Kernel:  AVX+FMA - 8 threads (8 CPU cores detected)

Reading alignment file example.phy ... Phylip format detected
Alignment most likely contains DNA/RNA sequences
Alignment has 17 sequences with 1998 columns, 1152 distinct patterns
1009 parsimony-informative, 303 singleton sites, 686 constant sites
[...]
Total number of iterations: 102
CPU time used for tree search: 34.761 sec (0h:0m:34s)
Wall-clock time used for tree search: 4.757 sec (0h:0m:4s)
Total CPU time used: 44.325 sec (0h:0m:44s)
Total wall-clock time used: 6.151 sec (0h:0m:6s)

Analysis results written to:
  IQ-TREE report:                example.phy.iqtree
  Maximum-likelihood tree:       example.phy.treefile
  Likelihood distances:          example.phy.mldist
  Screen log file:               example.phy.log

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. iqtree.sh). For example:



```

#!/bin/bash
module load iqtree
iqtree2 -s example.phy -nt $SLURM_CPUS_PER_TASK

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=8 --mem=6g iqtree.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. iqtree.swarm). For example:



```

iqtree2 -s example1.phy -nt $SLURM_CPUS_PER_TASK
iqtree2 -s example2.phy -nt $SLURM_CPUS_PER_TASK
iqtree2 -s example3.phy -nt $SLURM_CPUS_PER_TASK

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f iqtree.swarm -g 6 -t 8 --module iqtree
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module iqtree Loads the iqtree module for each subjob in the swarm 
 | |
 | |
 | |








