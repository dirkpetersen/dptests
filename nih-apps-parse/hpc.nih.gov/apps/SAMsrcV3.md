

document.querySelector('title').textContent = 'SAMsrcV3 on Biowulf';
SAMsrcV3 on Biowulf


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


Synthetic Aperture Magnetometry
-------------------------------



The SAMsrcV3 suite implements the latest advances in MEG source localization. Compared with the older CTF tools, there are many new features (some of them experimental). Filtering options, more flexible control over active, control, and noise covariance matrices, dual-state imaging, realistic head models, and several flavors of multipole source imaging are available.




Documentation
* [SAMsrcV3 Main Site](https://kurage.nimh.nih.gov/meglab/Meg/SAMsrcV3)


Important Notes
* Module Name: samsrcv3 (see [the modules page](/apps/modules.html) for more information)



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

[user@cn3144 ~]$ **module load samsrcv3**

[user@cn3144 ~]$ **sam\_cov -r my\_dataset -m settings.param**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. samsrc.sh). For example:



```

#!/bin/bash
set -e
module load samsrcv3
sam_cov -r my_dataset -m settings.param

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] samsrc.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. samsrc.swarm). For example:



```

sam_cov -r dataset01 -m settings.param
sam_cov -r dataset02 -m settings.param
sam_cov -r dataset03 -m settings.param
sam_cov -r dataset04 -m settings.param

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f samsrc.swarm [-g #] [-t #] --module samsrcv3
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module samsrcv3 Loads the SAMsrcV3 module for each subjob in the swarm 
 | |
 | |
 | |








