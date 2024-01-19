

document.querySelector('title').textContent = 'AFNI on Biowulf';
AFNI on Biowulf


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



AFNI (Analysis of Functional NeuroImages) is a set of C programs for
processing, analyzing, and displaying functional MRI (FMRI) data - a technique
for mapping human brain activity. AFNI is developed by the Scientific and
Statistical Computing Core, NIMH. 



Documentation
* [AFNI website](http://afni.nimh.nih.gov/afni/)* [AFNI documentation](https://afni.nimh.nih.gov/pub/dist/doc/htmldoc/index.html)


Important Notes
* Module Name: AFNI (see [the modules page](/apps/modules.html) for more information)
* AFNI
binaries on Biowulf are updated once a day at 4 am, so loading the 'current-openmp' version will
get you the latest binaries. Loading the module will give you the timestamp of the last update, e.g.

```

[user@biowulf ~]$ **module load afni/current-openmp**
[+] Loading AFNI current-openmp ...
AFNI/current-openmp last updated  **2017-11-01**

```
* The default version of AFNI on Biowulf contains some multi-threaded executables. AFNI jobs using multi-threaded executables can be submitted to more than the default
2 CPUs, as described below.
* A copy of AFNI is retained every 3 months, for users who want to complete a run using the same version of the executables. These 3-month snapshots are stored for a year.
* This application produces HTML reports. You can use [hpcdrive to view these reports](/docs/hpcdrive.html) on your local workstation.



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

[user@cn3144 ~]$ **module load afni**
[+] Loading AFNI current-openmp ...
AFNI/current-openmp last updated  2017-05-22

[user@cn3144 ~]$ **cd /data/user/afni/AFNI\_data2\_helix/**

[user@cn3144 AFNI_data2_helix]$ **time tcsh proc.subj ED**
(version 1.24, Jun 4, 2007)
++ 3dcopy: AFNI version=AFNI_2011_12_21_1014 (Jan  8 2015) [64-bit]
++ 3dTcat: AFNI version=AFNI_2011_12_21_1014 (Jan  8 2015) [64-bit]
++ 3dTcat: AFNI version=AFNI_2011_12_21_1014 (Jan  8 2015) [64-bit]
++ 3dTcat: AFNI version=AFNI_2011_12_21_1014 (Jan  8 2015) [64-bit]
++ 3dTcat: AFNI version=AFNI_2011_12_21_1014 (Jan  8 2015) [64-bit]

[....etc...]

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. AFNI.sh). For example:



```

#!/bin/bash
# this file is called AFNI.sh

module load AFNI
cd /data/$USER/AFNI_data2_helix
tcsh ./proc.subj ED

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] AFNI.sh
```

You would use --cpus-per-task for the multithreaded AFNI executables. --mem is used if you need more than the default memory.

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.

Typically, AFNI on Biowulf is used to simultaneously process a large number
of independent datasets via the [swarm](/apps/swarm.html)
utility.


A small sample set of data for 3 subjects can be copied from
/usr/local/apps/afni/AFNI\_data2\_helix.tar.gz,
courtesy of Rick Reynolds (NIMH). This tar file includes a script called
proc.subj to process data for a single subject.


You can unpack this dataset into your own area with

```

cd /data/$USER
tar xvzf /usr/local/apps/afni/AFNI_data2_helix.tar.gz

```

A swarm command file to process these 3 datasets would look like




```

#-------- this file is called swarm.cmd -----------------
cd /data/$USER/AFNI_data2_helix; tcsh ./proc.subj ED
cd /data/$USER/AFNI_data2_helix; tcsh ./proc.subj EE
cd /data/$USER/AFNI_data2_helix; tcsh ./proc.subj EF
#--------------------------------------------------------

```


This swarm command file would be submitted to the batch system with:



```

swarm -f swarm.cmd --module afni

```


If each individual process requires more than 4 GB of RAM, you can
specify the required memory with

```

swarm -g # -f swarm.cmd --module afni

```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g #  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t #  Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module afni Loads the AFNI module for each subjob in the swarm 
 | |
 | |
 | |










