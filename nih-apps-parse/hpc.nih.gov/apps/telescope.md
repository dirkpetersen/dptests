

document.querySelector('title').textContent = 'Telescope on Biowulf';
Telescope on Biowulf


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



 Single locus resolution of Transposable ELEment expression.
 Telescope estimates transposable element expression (retrotranscriptome) resolved to specific genomic locations. It directly addresses uncertainty in fragment assignment by reassigning ambiguously mapped fragments to the most probable source transcript as determined within a Bayesian statistical model. 



### References:


* Matthew L Bendall, Miguel de Mulder, Luis Pedro Iñiguez, Aarón Lecanda-Sánchez, Marcos Pérez-Losada, Mario A Ostrowski, Richard B Jones, Lubbertus Mulder, Gustavo Reyes-Terán, Keith A Crandall, Christopher E Ormsby, Douglas F. Nixon.
 **Telescope: Characterization of the retrotranscriptome by accurate estimation of transposable element expression.**
 bioRxiv 398172; doi: [10.1101/398172](https://doi.org/10.1101/398172)


Documentation
* [Telescope Main Site](https://github.com/mlbendall/telescope)


Important Notes
* Module Name: telescope (see [the modules page](/apps/modules.html) for more information)



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

[user@cn3144 ~]$ **module load telescope**
[+] Loading telescope, version 1.0.2... 
[user@cn3144 ~]$ **telescope assign *[samfile]* *[gtffile]***
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. telescope.sh). For example:



```

#!/bin/bash
set -e
module load telescope
telescope assign [samfile] [gtffile]

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] telescope.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. telescope.swarm). For example:



```

telescope assign [samfile1] [gtffile]
telescope assign [samfile2] [gtffile]
telescope assign [samfile3] [gtffile]
telescope assign [samfile4] [gtffile]

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f telescope.swarm [-g #] [-t #] --module telescope
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module telescope Loads the Telescope module for each subjob in the swarm 
 | |
 | |
 | |








