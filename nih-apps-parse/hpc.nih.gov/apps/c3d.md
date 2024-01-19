

document.querySelector('title').textContent = 'C3D on Biowulf';
C3D on Biowulf


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



C3D is a command-line tool for converting 3D images between common file formats. The tool also includes a growing list of commands for image manipulation, such as thresholding and resampling.



### Web site


* [Home page](https://sourceforge.net/projects/c3d/)


Documentation
* [C3D Documentation](https://sourceforge.net/p/c3d/git/ci/master/tree/doc/c3d.md)


Important Notes
* Module Name: c3d (see [the modules page](/apps/modules.html) for more information)



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

[user@cn3144 ~]$ **module load c3d**

[user@cn3144 ~]$ **c3d -h**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. c3d\_job.sh). For example:



```

#!/bin/bash
set -e
module load c3d
c3d minimal.nii output.img

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] c3d_job.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. c3d\_jobs.swarm). For example:



```

c3d minimal.nii -interpolation Cubic -resample 50x30x40vox -o output_0.img
c3d minimal.nii -interpolation Cubic -resample 25x15x20vox -o output_1.img
c3d minimal.nii -interpolation Cubic -resample 12x7x10vox -o output_2.img

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f c3d_jobs.swarm [-g #] [-t #] --module c3d
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module c3d Loads the cmtk module for each subjob in the swarm 
 | |
 | |
 | |








