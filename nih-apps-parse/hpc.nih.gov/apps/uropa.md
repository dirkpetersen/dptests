

document.querySelector('title').textContent = 'UROPA on Biowulf';
UROPA on Biowulf


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



UROPA is a command line based tool intended for genomic region annotation.



### References:


* Kondili M, Fust A, Preussner J, Kuenne C, Braun T, Looso M. [**UROPA: a tool for Universal RObust Peak Annotation.**](https://www.ncbi.nlm.nih.gov/pubmed/28572580) Scientific Reports. 2017;7:2593. doi:10.1038/s41598-017-02464-y


Documentation
* [UROPA Main Site](https://uropa-manual.readthedocs.io/index.html)


Important Notes
* Module Name: UROPA (see [the modules page](/apps/modules.html) for more information)



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

[user@cn1234 ~]$ **module load uropa**

[user@cn1234 ~]$ **uropa -h**
[+] Loading singularity  on cn3397 
[+] Loading uropa 2.0.3  ... 
[user@cn1234 ~]$ **uropa -i test.json -p output/test -s**
[user@cn1234 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. uropa.sh). For example:



```

#!/bin/bash
set -e
module load uropa
uropa -i test.json -p output/test_gem -s

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] uropa.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. uropa.swarm). For example:



```

uropa -i test1.json -p output/test1_gem -s
uropa -i test2.json -p output/test2_gem -s
uropa -i test3.json -p output/test3_gem -s

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f uropa.swarm [-g #] [-t #] --module uropa
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module uropa Loads the uropa module for each subjob in the swarm 
 | |
 | |
 | |








