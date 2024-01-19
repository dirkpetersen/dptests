

document.querySelector('title').textContent = 'cellphonedb on Biowulf';
cellphonedb on Biowulf


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



A publicly available repository of curated receptors, ligands and their interactions. Subunit architecture is included for both ligands and receptors, representing heteromeric complexes accurately.



### References:


* [Vento-Tormo, Roser, et al. "Single-cell reconstruction of the early maternal–fetal interface in humans." *Nature* 563.7731 (2018): 347.](https://www.nature.com/articles/s41586-018-0698-6)* [Efremova, Mirjana, et al. "CellPhoneDB v2. 0: Inferring cell-cell communication from combined expression of multi-subunit receptor-ligand complexes." *bioRxiv* (2019): 680926.](https://www.biorxiv.org/content/10.1101/680926v1)


Documentation
* [cellphonedb Main Site](https://www.cellphonedb.org/)
* [cellphonedb GitHub repo](https://github.com/ventolab/cellphonedb)


Important Notes
* Module Name: cellphonedb (see [the modules page](/apps/modules.html) for more information)
 * Environment variables set 
	+ CELLPHONEDB\_HOME* cellphonedb uses the python package "Click" that is very opinionated about locale. If you receive locale related errors try setting your LANG and LC\_ALL variables to something like en\_US.utf8



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf ~]$ **sinteractive --cpus-per-task=4 --mem=16G --gres=lscratch:10**
salloc.exe: Pending job allocation 41538656
salloc.exe: job 41538656 queued and waiting for resources
salloc.exe: job 41538656 has been allocated resources
salloc.exe: Granted job allocation 41538656
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3114 are ready for job
srun: error: x11: no local DISPLAY defined, skipping

[user@cn3114 ~]$ **cd /data/$USER**

[user@cn3114 user]$ **git clone https://github.com/ventolab/cellphonedb.git**
Cloning into 'cellphonedb'...
remote: Enumerating objects: 15, done.
remote: Counting objects: 100% (15/15), done.
remote: Compressing objects: 100% (14/14), done.
remote: Total 9381 (delta 1), reused 4 (delta 1), pack-reused 9366
Receiving objects: 100% (9381/9381), 61.54 MiB | 24.54 MiB/s, done.
Resolving deltas: 100% (7066/7066), done.
Checking out files: 100% (368/368), done.

[user@cn3114 user]$ **cd cellphonedb/in/example\_data/**

[user@cn3114 example_data]$ **module load cellphonedb**
[+] Loading cellphonedb  2.1.1  on cn3114
[+] Loading gcc  7.3.0  ...
[+] Loading GSL 2.4 for GCC 7.2.0 ...
[-] Unloading gcc  7.3.0  ...
[+] Loading gcc  7.3.0  ...
[+] Loading openmpi 3.0.2  for GCC 7.3.0
[+] Loading ImageMagick  7.0.8  on cn3114
[+] Loading HDF5  1.10.4
[+] Loading pandoc  2.1.1  on cn3114
[+] Loading R 3.6.0

[user@cn3114 example_data]$ **cellphonedb method statistical\_analysis test\_meta.txt test\_counts.txt**
[ ][APP][08/11/19-18:25:15][WARNING] Latest local available version is `v2.0.0`, using it
[ ][APP][08/11/19-18:25:15][WARNING] User selected downloaded database `v2.0.0` is available, using it
[ ][CORE][08/11/19-18:25:15][INFO] Initializing SqlAlchemy CellPhoneDB Core
[...snip...]

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. cellphonedb.sh). For example:



```

#!/bin/bash
set -e
module load cellphonedb
cd /data/$USER/cellphonedb/in/example_data/
cellphonedb method statistical_analysis test_meta.txt test_counts.txt

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] cellphonedb.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. cellphonedb.swarm). For example:



```

cellphonedb method analysis 1meta.txt 1counts.txt 
cellphonedb method analysis 2meta.txt 2counts.txt 
cellphonedb method analysis 3meta.txt 3counts.txt 
cellphonedb method analysis 4meta.txt 4counts.txt 

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f cellphonedb.swarm [-g #] [-t #] --module cellphonedb
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module cellphonedb Loads the cellphonedb module for each subjob in the swarm 
 | |
 | |
 | |








