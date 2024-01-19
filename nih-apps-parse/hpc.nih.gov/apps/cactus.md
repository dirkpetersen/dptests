

document.querySelector('title').textContent = 'Cactus on Biowulf';
Cactus on Biowulf


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


 Cactus is a reference-free whole-genome multiple alignment program.



### References:


* Armstrong J, Hickey G, Diekhans M, Fiddes IT, Novak AM, Deran A, Fang Q, Xie D, Feng S, Stiller J, Genereux D, Johnson J, Marinescu VD, Alföldi J, Harris RS, Lindblad-Toh K, Haussler D, Karlsson E, Jarvis ED, Zhang G, Paten B.
*[Progressive Cactus is a multiple-genome aligner for the thousand-genome era.](https://pubmed.ncbi.nlm.nih.gov/33177663/)*  Nature. 2020 Nov;587(7833):246-251.


Documentation
* Cactus Main Site: [GitHub](https://github.com/ComparativeGenomicsToolkit/cactus)


Important Notes
* Module Name: cactus (see [the modules page](/apps/modules.html) for more information)



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --cpus-per-task=4 --gres=lscratch:8**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn4224 are ready for job

[user@cn4224 ~]$ **module load cactus**
[+] Loading cactus  1.2.3  on cn4224 
[+] Loading singularity  3.7.0  on cn4224

[user@cn4224 ~]$ **cd /lscratch/${SLURM\_JOB\_ID}**
[user@cn4224 ~]$ **export SINGULARITYENV\_TMPDIR=/lscratch/${SLURM\_JOB\_ID}**
[user@cn4224 ~]$ **wget https://raw.githubusercontent.com/ComparativeGenomicsToolkit/cactus/master/examples/evolverMammals.txt**

[user@cn4224 ~]$ **cactus /lscratch/${SLURM\_JOB\_ID}/jobStore \
 /lscratch/${SLURM\_JOB\_ID}/evolverMammals.txt \
 /data/${USER}/evolverMammals.hal \
 --root mr \
 --binariesMode local**
[...]
Workflow Progress 100%||||||||||||||||||||||| 379/379 (0 failures) [19:08<00:00, 0.33 jobs/s]
[2020-12-23T13:47:39-0500] [MainThread] [I] [toil.common] Successfully deleted the job store: FileJobStore(/lscratch/46116226/jobStore)
[2020-12-23T13:47:39-0500] [MainThread] [I] [cactus.progressive.cactus_progressive] Cactus has finished after 1156.1930517529836 seconds

[user@cn4224 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226

[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. cactus.sh) similar to the following.



```

#! /bin/bash

set -e

module load cactus

cd /lscratch/${SLURM_JOB_ID}
export SINGULARITYENV_TMPDIR=/lscratch/${SLURM_JOB_ID}
cactus /lscratch/${SLURM_JOB_ID}/jobStore \
       /lscratch/${SLURM_JOB_ID}/evolverMammals.txt \
       /data/${USER}/evolverMammals.hal \
       --root mr \
       --binariesMode local

```

Submit these jobs using the Slurm [sbatch](/docs/userguide.html) command.


Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile for the first step of the pipeline (e.g. cactus.swarm). For example:



```

cd /lscratch/${SLURM_JOB_ID}; \
export SINGULARITYENV_TMPDIR=/lscratch/${SLURM_JOB_ID}; \
cactus /lscratch/${SLURM_JOB_ID}/jobStore \
       /lscratch/${SLURM_JOB_ID}/file1.txt \
       /data/${USER}/file1.hal \
       --root mr \
       --binariesMode local
cd /lscratch/${SLURM_JOB_ID}; \
export SINGULARITYENV_TMPDIR=/lscratch/${SLURM_JOB_ID}; \
cactus /lscratch/${SLURM_JOB_ID}/jobStore \
       /lscratch/${SLURM_JOB_ID}/file2.txt \
       /data/${USER}/file2.hal \
       --root mr \
       --binariesMode local

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f cactus.swarm [-g #] --module cactus
```

where


|  |  |  |  |
| --- | --- | --- | --- |
| -g #  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | --module cactus  Loads the cactus module for each subjob in the swarm
 | |
 | |








