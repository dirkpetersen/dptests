

document.querySelector('title').textContent = 'Mango on Biowulf';
Mango on Biowulf


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



Mango (Multi-image Analysis GUI) is a viewer for medical research images. Mango contains analysis tools and a GUI to navigate image volumes. Mango also comes with the following command-line utilities: mango (loads image), mango-applytransform (applies transform to image), mango-convert2avw (converts to AVW format), mango-convert2des (converts to DES format), mango-convert2nii (converts to NIFTI format), mango-imageinfo (prints image metadata summary), mango-makeroi (makes ROI based on image threshold), mango-resizer (resizes image), mango-script (runs mango script), mango-vols2series (concatenates volumes into series). 



### Web site


* [Home page](http://ric.uthscsa.edu/mango/index.html)


Documentation
* [Mango User Guide](http://ric.uthscsa.edu/mango/userguide.html)


Important Notes
* Module Name: mango (see [the modules page](/apps/modules.html) for more information)



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

[user@cn3144 ~]$ **module load mango**

[user@cn3144 ~]$ **mango-imageinfo sub-01\_T1w.nii.gz**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. mango\_job.sh). For example:



```

#!/bin/bash

export TMPDIR=/lscratch/$SLURM_JOB_ID

module load mango
mango-convert2avw sub-01_T1w.nii.gz

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] mango_job.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. mango\_jobs.swarm). For example:



```

mango-convert2avw sub-01_T1w.nii.gz
mango-convert2avw sub-02_T1w.nii.gz
mango-convert2avw sub-03_T1w.nii.gz
mango-convert2avw sub-04_T1w.nii.gz
mango-convert2avw sub-05_T1w.nii.gz

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f mango_jobs.swarm [-g #] [-t #] --module mango
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module mango Loads the mango module for each subjob in the swarm 
 | |
 | |
 | |








