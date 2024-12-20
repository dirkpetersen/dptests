

document.querySelector('title').textContent = 'Naccess on Biowulf';
Naccess on Biowulf


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



Briefly, the naccess program calculates the atomic accessible surface defined
by rolling a probe of given size around a van der Waals surface. This program
is an implimentation of the method of Lee and Richards (1971) J.Mol.Biol.55,
379-400. which does just that. The program is dimensioned for up to 20000
atoms, and allows the variation of the probe size and atomic radii by the user.



Documentation
* Naccess Main Site: <http://www.bioinf.manchester.ac.uk/naccess/>


Important Notes
* Module Name: naccess (see [the modules page](/apps/modules.html) for more information)
* Singlethreaded
* Reference data in /pdb/pdb/



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

[user@cn3144 ~]$ ml naccess
[+] Loading naccess 2.1.1 ...
[user@cn3144 ~]$ zcat /pdb/pdb/cr/pdb1crn.ent.gz > 1crn.pdb
[user@cn3144 ~]$ naccess 1crn.pdb
naccess: using defualt vdw.radii
naccess: using default STD FILE

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. naccess.sh). For example:



```

#!/bin/bash
module load naccess
naccess /path/to/blah.pdb

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] naccess.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. naccess.swarm). For example:



```

naccess 1pdb.pdb
naccess 2pdb.pdb
naccess 3pdb.pdb
naccess 4pdb.pdb

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f naccess.swarm [-g #] [-t #] --module naccess
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g #  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t #  Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module naccess  Loads the naccess module for each subjob in the swarm 
 | |
 | |
 | |








