

document.querySelector('title').textContent = 'Mega2 on Biowulf';
Mega2 on Biowulf


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



“Manipulation Environment for Genetic Analyses”
A data-handling program for facilitating
genetic linkage and association analyses

Documentation
* <https://watson.hgen.pitt.edu/docs/mega2_html/mega2.html>


Important Notes
* Module Name: mega2 (see [the modules page](/apps/modules.html) for more information)



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive --mem=5g**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load mega2**
[user@cn3144 ~]$ **mega2 --help**
Usage: mega2 [options] [batch-file-name] {arguments}
  acceptable options:
  DB         --DBfile 
 change the database name from dbmega2.db to .
 --DBdump
 dump the database and if --DBread is also present,
 then exec’s a new copy of Mega2 to process the database.
 --DBread
 read an existing database file and do an analysis.
 --DBcompress 
 set database compression level: 0 == off; 1 == gzip; (default == 1).
...
...
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. mega2.sh). For example:



```

#!/bin/bash
set -e
module load mega2
cd /data/$USER
mega2 command

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --mem=5g mega2.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. mega2.swarm). For example:



```

cd dir1; mega2 command
cd dir2; mega2 command
...
cd dir10; mega2 command

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f mega2.swarm -g 5 --module mega2
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module mega2 Loads the mega2 module for each subjob in the swarm 
 | |
 | |
 | |












