

document.querySelector('title').textContent = "TEMPLATE";
vContact2 on Biowulf


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



Description



### References:


* Blow J.
 [**A really amazing research paper.**](https://www.ncbi.nlm.nih.gov/pubmed/00000000)
*J Mol Biol. 2012 Jan 13;415(2):406-18.*
* Blow J., Doe J.
 [**A retread of another amazing research paper.**](https://www.ncbi.nlm.nih.gov/pubmed/00000000)
*J Struct Biol. 2012 Dec;180(3):519-30.*


Documentation
* [TEMPLATE Main Site](https://hpcwebdev.cit.nih.gov/)


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








