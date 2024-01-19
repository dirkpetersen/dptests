

document.querySelector('title').textContent = 'Frealign on Biowulf';
Frealign on Biowulf


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


 Frealign is a program for high-resolution refinement of 3D reconstructions from cryo-EM images of single particles.


### References:


* Grigorieff, N. 2007. **[FREALIGN: high-resolution refinement of single particle structures.](http://dx.doi.org/10.1016/j.jsb.2006.05.004)** *J Struct Biol.* 157:117–125.


Documentation
* Frealign Main Site: [Frealign Main Page (Grigorieff Lab)](http://grigoriefflab.janelia.org/frealign)


Important Notes
* Module Name: Frealign (see [the modules page](/apps/modules.html) for more information)
* Multithreaded/Distributed
* environment variables set 
	+ FREALIGN\_HOME -- installation directory for Frealign
	+ NCPUS -- number of cpus/threads available for the job


The command frealign\_v9\_mp.exe is multithreaded. The number of threads is set by the **$NCPUS** environment variable. This variable is set dynamically by the Frealign module using the **$SLURM\_CPUS\_ON\_NODE** variable set by the Slurm batch system. It can be overridden by setting **$NCPUS** prior to loading the Frealign module.
  

To see the value of **$NCPUS** prior to running, type '**module show Frealign**'.


Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.
You will need to create an mparameters file (for example, frealign\_run\_refine).
Sample session:



```

[user@biowulf]$ **sinteractive**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ module load Frealign
[user@cn3144 ~]$ frealign_template
[user@cn3144 ~]$ nano mparameters
[user@cn3144 ~]$ frealign_run_refine

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
### Single Node Job


There are two ways of running a batch job on Biowulf. If your commands involve subtasks that run very quickly (for
example reconstruction steps that last a few minutes each), it is **much more efficient to run on a single node, using the local
scratch disk**. This is done by editing the mparameters file. In this case, because the path to the local scratch disk
space is unknown prior to the job submission, the mparameters file contains a dummy tag "XXXXXXXX".



```

cluster_type         none
nprocessor_ref       16
nprocessor_rec       16
scratch_dir          /lscratch/XXXXXXXX
```

Create a batch input file (e.g. frealign.sh) which cds to the directory containing the mparameters file and
substitutes the actual path to local disk into the dummy tag. For example:



```
#!/bin/bash
module load Frealign
cd /path/to/directory/where/mparameters/lives/
sed -i "s#/lscratch/XXXXXXXX#/lscratch/$SLURM_JOB_ID#" mparameters
frealign_run_refine
```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command, allocating as many cpus as indicated
in the nprocessor\_ref and/or nprocessor\_rec values, as well as local scratch space:



```
sbatch --cpus-per-task=16 --gres=lscratch:50 frealign.sh
```

### Multiple Node Job


In cases where the subtasks of Frealign as expected to run for a significantly long time (more than 10 minutes), it may
be better to distribute the subtasks across multiple nodes. Again, edit the mparameters file to indicate 'slurm' as the
cluster type, and set nprocessor\_{ref,rec} to no more than 50. In this case, you must use a shared scratch space, so either
leave scratch\_dir blank to indicate the current working directory, or set the value to a specific directory. The value for
mp\_cpus can be safely set to 2, because each sbatch job will allocate a minimum of 2 cpus.



```

cluster_type         slurm
nprocessor_ref       50
nprocessor_rec       50
mp_cpus              2

```

Create a batch script (for example frealign.sh):



```
#!/bin/bash
module load Frealign
cd /path/to/directory/where/mparameters/lives/
frealign_v9_mp.exe
```

Then submit to the cluster without any special allocations:



```
sbatch frealign.sh
```

### Job Control via mparameters File


When the cluster\_type is set to 'slurm', there are two values in mparameters which can give 
control the allocations for batch submission: qsub\_string\_ref and qsub\_string\_rec. These can be used to select
specific partitions, allocations, etc. For example, to run subtasks on a different partition:



```
qsub_string_ref      "--partition=quick"
qsub_string_rec      "--partition=quick""
```

If the number of cpus used during reconstruction can be larger than 2, then both mp\_cpus and qsub\_string\_rec need to
complement each other. For example, to set then number of cpus to 8:



```

mp_cpus              8
qsub_string_rec      "--cpus-per-task=8"

```





