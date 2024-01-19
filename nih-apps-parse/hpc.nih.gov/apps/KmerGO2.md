

document.querySelector('title').textContent = 'KmerGO2 on Biowulf';
KmerGO2 on Biowulf


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



KmerGO KmerGO is a user-friendly tool to identify the group-specific sequences on two groups of high throughput sequencing datasets. A sequence that is present, or rich, in one group, but absent, or scarce, in another group is considered “group-specific” here. Furthermore, KmerGO can also be applied to capture trait-associated sequences for continuous-trait dataset.



### References:


* Wang Y, Chen Q, Deng C, Zheng Y and Sun F
  [**KmerGO: A Tool to Identify Group-Specific Sequences With k-mers.c**](https://www.frontiersin.org/articles/10.3389/fmicb.2020.02067/full)
*Front. Microbiol. 11:2067. doi: 10.3389/fmicb.2020.02067.*


Documentation
* [KmerGO2 Main Site](https://github.com/ChnMasterOG/KmerGO2)


Important Notes
This application requires a [graphical connection using NX](/docs/connect.html#nx)


* Module Name: TEMPLATE (see [the modules page](/apps/modules.html) for more information)
 * Multithreaded/singlethreaded/MPI...
 * This application produces HTML reports. You can use [hpcdrive](/docs/hpcdrive.html) to view these reports on your local workstation.
* Environment variables set 
	+ TEMPLATE\_HOME* Example files in ???* Reference data in /fdb/TEMPLATE/



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

[user@cn3144 ~]$ **module load TEMPLATE**

[user@cn3144 ~]$ 

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. TEMPLATE.sh). For example:



```

#!/bin/bash
set -e
module load TEMPLATE
TEMPLATE < TEMPLATE.in > TEMPLATE.out

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] TEMPLATE.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. TEMPLATE.swarm). For example:



```

TEMPLATE < TEMPLATE.in > TEMPLATE.out
TEMPLATE < TEMPLATE.in > TEMPLATE.out
TEMPLATE < TEMPLATE.in > TEMPLATE.out
TEMPLATE < TEMPLATE.in > TEMPLATE.out

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f TEMPLATE.swarm [-g #] [-t #] --module TEMPLATE
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module TEMPLATE Loads the TEMPLATE module for each subjob in the swarm 
 | |
 | |
 | |








