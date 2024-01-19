

document.querySelector('title').textContent = 'Arriba on Biowulf';
Arriba on Biowulf


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



Arriba is a command-line tool for detecting gene fusions in RNA-Seq data. It
can also detect other clinically-relevant structural variations such as exon duplications or truncations of genes (i.e., breakpoints in introns and intergenic regions).



### References:


A dedicated publication about Arriba has not been released yet. Until then, please refer to Arriba in your methods section as follows (or similar):

 We used Arriba (https://github.com/suhrig/arriba/) to detect gene fusions from RNA-Seq data.


Documentation
* [Arriba Main Site](https://github.com/suhrig/arriba/)


Important Notes
* Module Name: arriba (see [the modules page](/apps/modules.html) for more information)
 * Arriba relies on the STAR genome aligner for much of its heavy lifting. Both can be run multi-threaded by setting the number of threads on the command-line to a value greater than 1.
 * The version of STAR encapsulated in arriba is 2.5.3a.
* Reference data in /fdb/arriba/references* Arriba is run by calling the run\_arriba.sh script, which does not use any flags to identify which argument is which. All arguments must be provided for each run, in the order specified by the script. For more details, see the information printed when run\_arriba.sh is called without any arguments.



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

[user@cn3144 ~]$ **module load arriba**
[+] Loading arriba  1.2.0  on cn3113
[+] Loading singularity  3.5.3  on cn3113
[user@cn3144 ~]$ **run\_arriba.sh**
Usage: run_arriba.sh STAR_genomeDir/ annotation.gtf assembly.fa blacklist.tsv read1.fastq.gz read2.fastq.gz threads
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. Arriba.sh). For example:



```

#!/bin/bash
set -e
module load arriba 
run_arriba.sh  > arriba.out

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] Arriba.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. Arriba.swarm). For example:



```

run_arriba.sh  > Arriba_1.out
run_arriba.sh  > Arriba_2.out
run_arriba.sh  > Arriba_3.out
run_arriba.sh  > Arriba_4.out

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f Arriba.swarm [-g #] [-t #] --module arriba
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module arriba Loads the arriba module for each subjob in the swarm 
 | |
 | |
 | |








