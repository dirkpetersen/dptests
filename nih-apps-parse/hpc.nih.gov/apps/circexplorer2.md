

document.querySelector('title').textContent = 'Circexplorer2 on Biowulf &amp; Helix';
Circexplorer2 on Biowulf & Helix



|  |
| --- |
| 
Quick Links
[Interactive job on Biowulf](#int)
[Batch job on Biowulf](#serial)
[Swarm of jobs](#swarm)
 |



Description
CIRCexplorer2 is a combined strategy to identify junction reads from back 
 spliced exons and intron lariats.


There may be multiple versions available on our systems. An easy way of 
 selecting the version is to use [modules](/apps/modules.html). 
 To see the modules available, type



```

module avail circexplorer2 

```

To select a module use



```

module load circexplorer2/[version]

```

where `[version]` is the version of choice.



Environment variables set
* `$PATH`
* `database files in /fdb/circexplorers2`


Documentation
<https://github.com/YangLab/CIRCexplorer>




Circexplorers2 job on Biowulf
Allocate an interactive session with [sinteractive](/docs/userguide.html#int) 
 and use as described below



```

biowulf$ **sinteractive --mem=5g**
salloc.exe: Pending job allocation 38978697
[...snip...]
salloc.exe: Nodes cn2273 are ready for job
node$ **module load circexplorer2**
[+] Loading circexplorer2
node$ **CIRCexplorer2.py -f accepted\_hits.bam -g hg19.fa -r ref.txt**
[...snip...]
node$ **exit**
biowulf$

```

 




Batch job on Biowulf
Create a batch script similar to the following example:



```

#! /bin/bash
# this file is circexplorer2.batch

module load circexplorer2 || exit 1
cd /data/$USER
CIRCexplorer2.py -f accepted_hits.bam -g hg19.fa -r ref.txt
```

Submit to the queue with [sbatch](/docs/userguide.html):



```

biowulf$ **sbatch circexplorer2.batch**

```

 



Swarm of Jobs on Biowulf
Create a swarmfile (e.g. script.swarm). For example:



```

# this file is called script.swarm
cd dir1;circexplorer2 command 1;circexplorer2 command 2
cd dir2;circexplorer2 command 1;circexplorer2 command 2
cd dir3;circexplorer2 command 1;circexplorer2 command 2
[...]

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f script.swarm --module circexplorer2
```

For more information regarding swarm: <https://hpc.nih.gov/apps/swarm.html#usage>




