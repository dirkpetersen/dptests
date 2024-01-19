

document.querySelector('title').textContent = 'shrimp on Biowulf';
shrimp on Biowulf


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



SHRiMP is a software package for aligning genomic reads against a target genome. It was primarily developed with the multitudinous short reads of next generation sequencing machines in mind, as well as Applied Biosystem's colourspace genomic representation. 



As of 2014, Shrimp is no longer being developed. The code is still available on Biowulf, but users are advised to migrate their
projects to use other software.

### References:


* SHRiMP was originally designed and written by
[Michael Brudno](http://www.cs.toronto.edu/~brudno) and Stephen M. Rumble,
with considerable input and testing by the [SidowLab](http://mendel.stanford.edu/SidowLab).
Since then, Adrian Dalca, Marc Fiume and Vladimir Yanovsky
have made considerable contributions to probability calculations and 2-pass SMS mapping algorithms.
The original SHRiMP publication can be found [here](http://dx.doi.org/10.1371/journal.pcbi.1000386).


Documentation
* [Shrimp documentation](http://compbio.cs.toronto.edu/shrimp/)


Important Notes
* Module Name: shrimp (see [the modules page](/apps/modules.html) for more information)
* Multithreaded/singlethreaded/MPI...
* Environment variables set 
	+ shrimp\_HOME* Example files in $SHRIMP\_DATA* Reference data in /fdb/shrimp/



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

[user@cn3144 ~]$ **module load shrimp**

[user@cn3144 ~]$ **cd /data/$USER**

[user@cn3144 ~]$ **cp $SHRIMP\_DATA/example/\* .** 

[user@cn3144 ~]$ **gmapper-cs test\_S1\_F3.csfasta ch11\_12\_validated.fasta -N 8 -o 5 -h 80% >map.out 2>map.log**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. shrimp.sh). For example:



```

#!/bin/bash
set -e
module load shrimp
cd /data/$USER/example
gmapper-cs test_S1_F3.csfasta ch11_12_validated.fasta -N 8 -o 5 -h 80% >map.out 2>map.log

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] shrimp.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. shrimp.swarm). For example:



```

gmapper-cs  file1.csfasta ch11_12_validated.fasta -N 8 -o 5 -h 80% >file1.out 2>&1
gmapper-cs  file2.csfasta ch11_12_validated.fasta -N 8 -o 5 -h 80% >file2.out 2>&1
gmapper-cs  file3.csfasta ch11_12_validated.fasta -N 8 -o 5 -h 80% >file3.out 2>&1
[...]

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f shrimp.swarm [-g #] [-t #] --module shrimp
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module shrimp Loads the shrimp module for each subjob in the swarm 
 | |
 | |
 | |










