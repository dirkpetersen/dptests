

document.querySelector('title').textContent = 'minimac on Biowulf';
minimac on Biowulf


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



Minimac is a low memory, computationally efficient implementation of the MaCH algorithm for genotype imputation. 



Documentation
* [Minimac4 website](https://genome.sph.umich.edu/wiki/Minimac4)


Important Notes
* Module Name: minimac (see [the modules page](/apps/modules.html) for more information)
* minimac-omp is multi-threaded. Plain 'minimac' is single-threaded.



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

[user@cn3144 ~]$ **module load minimac**
[+] Loading minimac minimac4  ...

[user@cn3144 ~]$  **minimac4 --refHaps refPanel.vcf \ 
 --haps targetStudy.vcf \
 --prefix testRun**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. minimac.sh). For example:



```

#!/bin/bash
set -e
module load minimac/minimac4
minimac-omp --cpus $SLURM_CPUS_PER_TASK --refHaps refPanel.vcf \ 
                --haps targetStudy.vcf \
                --prefix testRun

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=4 [--mem=#] minimac.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. minimac.swarm). For example:



```

minimac-omp --cpus $SLURM_CPUS_PER_TASK -refHaps ref1.vcf --haps target1.vcf --prefix test1
minimac-omp --cpus $SLURM_CPUS_PER_TASK -refHaps ref2.vcf --haps target2.vcf --prefix test2
minimac-omp --cpus $SLURM_CPUS_PER_TASK -refHaps ref3.vcf --haps target3.vcf --prefix test3

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f minimac.swarm [-g #] [-t 4] --module minimac
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file). 4 or 8 threads max is recommended.
 | --module minimac Loads the minimac module for each subjob in the swarm 
 | |
 | |
 | |








