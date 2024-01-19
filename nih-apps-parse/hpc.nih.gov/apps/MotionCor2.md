

document.querySelector('title').textContent = 'MotionCor2 on Biowulf';
MotionCor2 on Biowulf


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



MotionCor2 is a multi-GPU accelerated program which corrects anisotropic image motion at the single pixel level.
Iterative, patch-based motion detection is combined with spatial and temporal constraints and dose weighting.
MotionCor2 works on a wide range of data sets, including those very close to focus or with very short integration times, obviating the need for particle polishing. The application significantly improves Thon ring quality and 3D reconstruction resolution.



### References:


* Zheng SQ, Palovcak E, Armache JP, Verba KA, Cheng Y, Agard DA.
 [**MotionCor2: anisotropic correction of beam-induced motion for improved cryo-electron microscopy.**](https://www.ncbi.nlm.nih.gov/pubmed/28250466)
*Nat Methods. 2017 Apr;14(4):331-332.*


Documentation
* [MotionCor2 Main Site](http://msg.ucsf.edu/em/software/motioncor2.html)
* [MotionCor2 User Manual (version 1.1.0)](RELION/MotionCor2-UserManual-05-03-2018.pdf)
* [MotionCor2 User Manual (version 1.3.0)](RELION/MotionCor2-UserManual-10-22-2019.pdf)
* [MotionCor2 User Manual (version 1.5.0)](RELION/MotionCor2-UserManual-05-31-2022.pdf)
* Type MotionCor2 --help


Important Notes
* Module Name: MotionCor2 (see [the modules page](/apps/modules.html) for more information)
* Singlethreaded, GPU-accelerated
* Environment variables set 
	+ MOTIONCOR2\_HOME
	+ RELION\_MOTIONCORR\_EXECUTABLE
	+ RELION\_MOTIONCOR2\_EXECUTABLE* Example files in /fdb/app\_testdata/cryoEM



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --constraint=**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load MotionCor2**
[user@cn3144 ~]$ MotionCor2 -InMrc mymovie.mrcs -OutMrc mymicrograph.mrc -LogFile output.log -Bft 150 -PixSize 3.49 -OutStack 1 -Patch 1 1 -Trunc 1 -Gpu 0

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. MotionCor2.sh). For example:



```

#!/bin/bash
module load MotionCor2
MotionCor2 -InMrc /path/to/raw/data/12345.mrcs -OutMrc Micrographs/12345.mrc -LogFile output.log -Bft 150 -PixSize 3.49 -OutStack 1 -Patch 1 1 -Trunc 1 -Gpu 0

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] --partition=gpu --gres=gpu:p100:1 MotionCor2.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. MotionCor2.swarm). For example:



```

MotionCor2 -InMrc /path/to/raw/data/0001.mrcs -OutMrc Micrographs/0001.mrc -LogFile output.log -Bft 150 -PixSize 3.49 -OutStack 1 -Patch 1 1 -Trunc 1 -Gpu 0
MotionCor2 -InMrc /path/to/raw/data/0002.mrcs -OutMrc Micrographs/0002.mrc -LogFile output.log -Bft 150 -PixSize 3.49 -OutStack 1 -Patch 1 1 -Trunc 1 -Gpu 0
MotionCor2 -InMrc /path/to/raw/data/0003.mrcs -OutMrc Micrographs/0003.mrc -LogFile output.log -Bft 150 -PixSize 3.49 -OutStack 1 -Patch 1 1 -Trunc 1 -Gpu 0
MotionCor2 -InMrc /path/to/raw/data/0004.mrcs -OutMrc Micrographs/0004.mrc -LogFile output.log -Bft 150 -PixSize 3.49 -OutStack 1 -Patch 1 1 -Trunc 1 -Gpu 0

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f MotionCor2.swarm [-g #] [-t #] --partition=gpu --gres=gpu:p100:1 --module MotionCor2
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module MotionCor2 Loads the MotionCor2 module for each subjob in the swarm 
 | |
 | |
 | |








