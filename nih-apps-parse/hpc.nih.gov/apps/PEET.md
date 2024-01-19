

document.querySelector('title').textContent = 'PEET on Biowulf';
PEET on Biowulf


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



PEET is a set of programs, separate from but typically used in conjunction with IMOD, to align and average particles / subvolumes extracted from 3D volumes.



### References:


* PEET was developed at the University of Colorado, Boulder. PEET
was first applied for averaging axonemes and its methods were described in 
[Nicastro et al 2006](http://bio3d.colorado.edu/PEET/PEETBibliography.html#Nicastro2006).
Many new features have since been added, including clustering / classification
as described in 
[Heumann et al 2011](http://bio3d.colorado.edu/PEET/PEETBibliography.html#Heumann2011).
Please cite these papers if you use PEET in your research. /li>


Documentation
* [PEET Manual](http://bio3d.colorado.edu/PEET/PEETmanual.html)


Important Notes
* Module Name: PEET (see [the modules page](/apps/modules.html) for more information)
* Used in conjunction with [IMOD](IMOD.html)* Environment variables set 
	+ PARTICLE\_DIR



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

[user@cn3144 ~]$ **module load PEET**

[user@cn3144 ~]$  **averageAll**
Starting averageAll_mce ...
This is PEET Version 1.12.0 5-March-2018.
Copyright 2000-2018 The Regents of the University of Colorado.
MATLAB. Copyright 1984-2018 The Mathworks, Inc.
[...]

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. PEET.sh). For example:



```

#!/bin/bash
set -e
module load PEET
PEET < PEET.in > PEET.out

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] PEET.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. PEET.swarm). For example:



```

PEET < PEET.in > PEET.out
PEET < PEET.in > PEET.out
PEET < PEET.in > PEET.out
PEET < PEET.in > PEET.out

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f PEET.swarm [-g #] [-t #] --module PEET
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module PEET Loads the PEET module for each subjob in the swarm 
 | |
 | |
 | |








