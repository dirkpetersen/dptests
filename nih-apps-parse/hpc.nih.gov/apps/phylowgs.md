

document.querySelector('title').textContent = 'Phylowgs on Biowulf';
Phylowgs on Biowulf


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



A new method, PhyloWGS, can be applied to whole-genome sequencing data from one or more tumor samples to reconstruct complete genotypes of these subpopulations based on variant allele frequencies (VAFs) of point mutations and population frequencies of structural variations. 



### References:


* Amit G Deshwar, Shankar Vembu, Christina K Yung, Gun Ho Jang, Lincoln Stein and Quaid Morris. PhyloWGS: Reconstructing subclonal composition and evolution from whole-genome sequencing of tumors. Genome Biology 2015 16:35


Documentation
* [Phylowgs Main Site](https://github.com/morrislab/phylowgs)


Important Notes
* Module Name: phylowgs (see [the modules page](/apps/modules.html) for more information)
* Singlethreaded app
* Phylowgs results can be viewed in local browser using [SSH Tunneling](https://hpc.nih.gov/docs/tunneling/)* Example files in /usr/local/apps/phylowgs/TEST\_DATA



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --gres=lscratch:20** 
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ mkdir -p /data/$USER/phylowgs && cd /data/$USER/phylowgs

[user@cn3144 ~]$ module load phylowgs

[user@cn3144 ~]$ cp /usr/local/apps/phylowgs/TEST_DATA/* .

[user@cn3144 ~]$ evolve.py ssm_data.txt cnv_data.txt -t /lscratch/$SLURM_JOBID
[2018-02-16 09:19:38.999007] Starting MCMC run...
[2018-02-16 09:19:39.004639] -1000
[2018-02-16 09:19:39.301332] Shrinking MH proposals. Now 200.000000
[2018-02-16 09:19:39.311407] -999
[2018-02-16 09:19:39.646055] Shrinking MH proposals. Now 400.000000
[2018-02-16 09:19:39.652818] -998
...

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. phylowgs.sh). For example:



```

#!/bin/bash
set -e
module load phylowgs
evolve.py ssm_data.txt cnv_data.txt -t /lscratch/$SLURM_JOBID

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --gres=lscratch:20 phylowgs.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. phylowgs.swarm). For example:



```

evolve.py dir1/ssm_data1.txt dir1/cnv_data1.txt -t /lscratch/$SLURM_JOBID
evolve.py dir2/ssm_data2.txt dir2/cnv_data2.txt -t /lscratch/$SLURM_JOBID
evolve.py dir3/ssm_data3.txt dir3/cnv_data3.txt -t /lscratch/$SLURM_JOBID

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f phylowgs.swarm --gres=lscratch:20 -g 5 -t 2 --module phylowgs
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module phylowgs Loads the phylowgs module for each subjob in the swarm 
 | |
 | |
 | |








