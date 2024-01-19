

document.querySelector('title').textContent = 'pigz on Biowulf';
pigz on Biowulf


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


Pigz is a multi-core implementation of gzip.



Documentation
* pigz Main Site: [GitHub](https://github.com/madler/pigz)


Important Notes
* Module Name: pigz (see [the modules page](/apps/modules.html) for more information)
* Environment variables set: *PIGZ\_HOME*.
 * Ensure that the number of CPUs allocated is equal or greater than the threads used in the *pigz* command. 
 * Note that to be able to use the environment variable *SLURM\_MEM\_PER\_NODE* below, you will need to explicitly allocate memory resources.



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --cpus-per-task=4 --mem=4g**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn4224 are ready for job

[user@cn4224 ~]$ **module load pigz**
[+] Loading pigz 2.7 on cn4224 

[user@cn4224 ~]$ **pigz --best -p 4 -v longshot\_output.vcf**
longshot_output.vcf to longshot_output.vcf.gz 

[user@cn4224 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226

[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. pigz.sh) similar to the following.



```

#! /bin/bash

set -e

module load pigz

pigz --best -p 4 longshot_output.vcf

```

Submit these jobs using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=4 --mem=4g pigz.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile for the first step of the pipeline (e.g. pigz.swarm). For example:



```

pigz --best -p 4 longshot_output_1.vcf
pigz --best -p 4 longshot_output_2.vcf
pigz --best -p 4 longshot_output_3.vcf

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f pigz.swarm -g 4 -t 4 --module pigz
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g #  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t #  Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module pigz  Loads the pigz module for each subjob in the swarm

 | |
 | |
 | |








