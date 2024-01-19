

document.querySelector('title').textContent = 'Subtom on Biowulf';
Subtom on Biowulf


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


 Subtom is a pipeline for subvolume alignment and averaging of electron cryo-tomography data. It uses the TOM toolbox (aquisition and analysis of electron tomography). 


Documentation
* Subtom Main Site: [Documentation](https://subtom.readthedocs.io/en/latest/)
* Subtom Download Site: [GitHub](https://github.com/DustinMorado/subTOM)


Important Notes
* Module Name: subtom (see [the modules page](/apps/modules.html) for more information)
* Scripts can found in `${SUBTOM_SCRIPTS}`



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
salloc.exe: Nodes cn4224 are ready for job

[user@cn4224 ~]$ **cd /data/${USER}**
[user@cn4224 ~]$ **mkdir SUBTOM\_SCRIPTS**
[user@cn4224 ~]$ **module load subtom**
[user@cn4224 ~]$ **cd SUBTOM\_SCRIPTS**
[user@cn4224 ~]$ **cp ${SUBTOM\_SCRIPTS}/subtom\_average.sh .**

```

At this step you need to (1) modify **subtom\_average.sh** to include the path to your datasets and (2) assign **run\_local="1"**



```

[user@cn4224 ~]$ **chmod u+x subtom\_average.sh**
[user@cn4224 ~]$ **./subtom\_average.sh**
Options sourced!

STARTING Averaging - Iteration: 1

STATUS Update: Averaging - Iteration: 1

	0 parallel sums out of 1

biowulf.nih.gov

ERROR Update: Averaging - Iteration: 1

LOG Update: Averaging - Iteration: 1

Summing Batch 1: [####################################### ] - 97% - 840 particles 0:00:50 | 0:00:52
Summing Batch 1: [####################################### ] - 97% - 844 particles 0:00:50 | 0:00:52
Summing Batch 1: [####################################### ] - 97% - 848 particles 0:00:50 | 0:00:52
Summing Batch 1: [####################################### ] - 98% - 852 particles 0:00:51 | 0:00:52
Summing Batch 1: [####################################### ] - 98% - 856 particles 0:00:51 | 0:00:52
Summing Batch 1: [########################################] - 99% - 860 particles 0:00:51 | 0:00:52
Summing Batch 1: [########################################] - 99% - 864 particles 0:00:51 | 0:00:52
Summing Batch 1: [########################################] - 100% - 868 particles 0:00:52 | 0:00:52
Summing Batch 1: [########################################] - 100% - 870 particles 0:00:52 | 0:00:52

STATUS Update: Averaging - Iteration: 1

	1 parallel sums out of 1

FINISHED Averaging - Iteration: 1

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. subtom.sh) similar to the following.



```

#! /bin/bash

set -e

module load subtom

/data/$USER/SUBTOM_SCRIPTS/subtom_average.sh

```

Submit these jobs using the Slurm [sbatch](/docs/userguide.html) command.


Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile for the first step of the pipeline (e.g. subtom.swarm). For example:



```

/data/$USER/SUBTOM_SCRIPTS/subtom_average1.sh
/data/$USER/SUBTOM_SCRIPTS/subtom_average2.sh
/data/$USER/SUBTOM_SCRIPTS/subtom_average3.sh
/data/$USER/SUBTOM_SCRIPTS/subtom_average4.sh
/data/$USER/SUBTOM_SCRIPTS/subtom_average5.sh

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f subtom.swarm [-g #] --module subtom
```

where


|  |  |  |  |
| --- | --- | --- | --- |
| -g #  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | --module subtom  Loads the subtom module for each subjob in the swarm
 | |
 | |








