

document.querySelector('title').textContent = 'CMTK on Biowulf';
CMTK on Biowulf


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



CMTK is a set of command-line tools for computational morphometry of biomedical images.
CMTK contains tools for (1) registration (affine and nonrigid, single and multichannel,
pairwise and groupwise); (2) image correction (MR bias field estimation, interleaved
image artifact correction, EPI unwarping); (3) processing (filters, combination of 
segmentations, shape-based averaging); and (4) statistics (t-tests and the general 
linear model).



### Web site


* [Home page](https://www.nitrc.org/projects/cmtk/)


Documentation
* [CMTK Documents](https://www.nitrc.org/docman/?group_id=212)


Important Notes
* Module Name: cmtk (see [the modules page](/apps/modules.html) for more information)



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

[user@cn3144 ~]$ **module load cmtk**

[user@cn3144 ~]$ 

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. cmtk\_job.sh). For example:



```

#!/bin/bash
set -e
module load cmtk
cmtk registrationx --dofs 12 --min-stepsize 1 -o 2_rigid.xform minimal.nii minimal2.nii

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] cmtk_job.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. cmtk\_jobs.swarm). For example:



```

cmtk registrationx --dofs 12 --min-stepsize 1 -o 2_rigid.xform minimal.nii minimal2.nii
cmtk registrationx --dofs 12 --min-stepsize 1 -o 2_rigid.xform minimal.nii minimal2.nii
cmtk registrationx --dofs 12 --min-stepsize 1 -o 2_rigid.xform minimal.nii minimal2.nii
cmtk registrationx --dofs 12 --min-stepsize 1 -o 2_rigid.xform minimal.nii minimal2.nii

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f cmtk_jobs.swarm [-g #] [-t #] --module cmtk
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module cmtk Loads the cmtk module for each subjob in the swarm 
 | |
 | |
 | |








