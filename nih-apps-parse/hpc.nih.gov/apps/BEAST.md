

document.querySelector('title').textContent = 'BEAST on Biowulf';
BEAST on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int)
[Batch job](#sbatch)
[Running on GPUs](#gpu)
[Swarm of jobs](#swarm)
 |



[BEAST](http://beast2.org/) (Bayesian Evolutionary Analysis
 Sampling Trees) is a cross-platform program for Bayesian MCMC analysis of
 molecular sequences. It is entirely orientated towards rooted, time-measured
 phylogenies inferred using strict or relaxed molecular clock models. It can be
 used as a method of reconstructing phylogenies but is also a framework for
 testing evolutionary hypotheses without conditioning on a single tree topology.
 BEAST uses MCMC to average over tree space, so that each tree is weighted
 proportional to its posterior probability. It includes a simple to use
 user-interface program for setting up standard analyses and a suit of programs
 for analysing the results.



### References:


* If you use this application in your research, the recommended citation is [Drummond AJ, Suchard MA, Xie D & Rambaut A (2012) Bayesian phylogenetics with BEAUti and the BEAST 1.7 Molecular Biology And Evolution 29: 1969-1973](http://mbe.oxfordjournals.org/content/29/8/1969)


Documentation
* Type **beast -help** at the prompt
* [BEAST1](http://beast.community/)
* [BEAST2](http://beast2.org/)
* [BEAST2 vs 1 Comparison](https://www.beast2.org/features/)


Important Notes
* Module Name: BEAST (see [the modules page](/apps/modules.html) for more information)
* Multithreaded and GPU capable
* environment variables set 
	+ BEAST\_LIB
	+ LD\_LIBRARY\_PATH
	+ BEAST\_EXAMPLES* Example files in $BEAST\_EXAMPLES



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

[user@cn3144 ~]$ mkdir /data/$USER/BEAST; cd /data/$USER/BEAST
[user@cn3144 BEAST]$ module load BEAST
[user@cn3144 BEAST]$ cp $BEAST_EXAMPLES/benchmark1.xml .
[user@cn3144 BEAST]$ beast benchmark1.xml

[user@cn3144 BEAST]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. BEAST.sh). For example:



```

#!/bin/bash
module load BEAST
beast -threads $SLURM_CPUS_PER_TASK myinput.xml

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] BEAST.sh
```

Running on GPUs
Under some circumstances, GPUs can greatly accelerate computation. Include -beagle\_GPU -beagle\_order 2,1,0 with the beast commandline:



```
beast -beagle_GPU -beagle_order 2,1,0 myinput.xml
```

and be sure to allocate GPUs:



```
sbatch --constraint=gpuk80 --gres=gpu:k80:2 sbatch.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. BEAST.swarm). For example:



```

beast -threads $SLURM_CPUS_PER_TASK -prefix input1 input1.xml
beast -threads $SLURM_CPUS_PER_TASK -prefix input2 input2.xml
beast -threads $SLURM_CPUS_PER_TASK -prefix input3 input3.xml
beast -threads $SLURM_CPUS_PER_TASK -prefix input4 input4.xml

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f BEAST.swarm [-g #] [-t #] --module BEAST
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module BEAST Loads the BEAST module for each subjob in the swarm 
 | |
 | |
 | |








