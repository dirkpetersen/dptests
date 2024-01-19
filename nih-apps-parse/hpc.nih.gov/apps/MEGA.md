

document.querySelector('title').textContent = 'MEGA on Biowulf';
MEGA on Biowulf


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



MEGA, or Molecular Evolutionary Genetic Analysis, is a suite of tools that help explore and analyze DNA and protein sequences for evolutionary and population genetics. These include methods for multiple sequence alignment, tree inference, time trees, models of evolution, diversity metrics, selection metrics, ancestral state reconstruction, and evolutionary rates.




The application can be run both as a GUI (command mega) or a command line interface (command megacc). To run the GUI, we recommend connecting to Biowulf using [NoMachine](docs/nx.html). The GUI can be used to run analyses and additionally create analysis options files (.mao files) that can be passed to the command-line application. This can be helpful for running large analyses with sbatch or swarm.



### Reference:


* Koichiro Tamura, Glen Stecher, and Sudhir Kumar
 [**MEGA11: Molecular Evolutionary Genetics Analysis version 11**](https://doi.org/10.1093/molbev/msab120) *Molecular Biology and Evolution 38:3022-3027 (2021)*


Documentation
* [MEGA Main Site](https://www.megasoftware.net/)
* [MEGA Online Manual](https://www.megasoftware.net/web_help_11/index.htm)


Important Notes
This application requires a [graphical connection using NX](/docs/connect.html#nx) to run the GUI version


* Module Name: MEGA (see [the modules page](/apps/modules.html) for more information)
* Singlthreaded. Multithreaded for ML phylogeny method.
* Run mega command to start the graphical user interface; megacc for the CLI
* Example files in $MEGA\_EXAMPLES. When you start the GUI, example files will automatically copied to "/home/$USER/MEGA X"
* When running the GUI, you can safely ignore the following message:   

libGL error: MESA-LOADER: failed to open swrast: /usr/lib64/dri/swrast\_dri.so: cannot open shared object file: No such file or directory (search paths /usr/lib64/dri)
 libGL error: failed to load driver: swrast


The Blast and Genbank features in MEGA GUI will not work on Biowulf since compute nodes are behind a firewall and MEGA is not able to use proxies. We recommend installing the GUI on your desktop to access these features. You can save a MEGA data session file and transfer it Biowulf to reopen in MEGA in an sinteractive session



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

[user@cn3144 ~]$ **module load MEGA**

[user@cn3144 ~]$ **mega**
![MEGA GUI](/images/mega-gui.png)


```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. MEGA.sh). For example:



```

#!/bin/bash
set -e
module load MEGA
megacc -a analysis.mao -d data.meg

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] MEGA.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. MEGA.swarm). For example:



```

megacc -a analysis1.mao -d data1.meg
megacc -a analysis2.mao -d data1.meg
megacc -a analysis1.mao -d data2.meg
megacc -a analysis2.mao -d data2.meg

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f MEGA.swarm [-g #] [-t #] --module MEGA
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module MEGA Loads the MEGA module for each subjob in the swarm 
 | |
 | |
 | |








