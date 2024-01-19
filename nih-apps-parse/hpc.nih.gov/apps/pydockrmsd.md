

document.querySelector('title').textContent = 'pyDockRMSD on Biowulf';
pyDockRMSD on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
 |



Description



### References:


* Eric W. Bell & Yang Zhang 
 [**DockRMSD: an open-source tool for atom mapping and RMSD calculation of symmetric molecules through graph isomorphism.**](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-019-0362-7)
*Journal of Cheminformatics volume 11, Article number: 40 (2019)*


Documentation
* [pyDockRMSD Main Site](https://github.com/neudinger/pyDockRMSD)


Important Notes
* Module Name: pydockrmsd (see [the modules page](/apps/modules.html) for more information)
 * singlethreaded
 * Example files in $PYDOCKRMSD\_EXAMPLES



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

[user@cn3144 ~]$ **module load pydockrmsd**

[user@cn3144 ~]$ 

[user@cn3144 ~]$ **mkdir test && cd test**
[user@cn3144 ~]$ **cp -Rp ../github/pyDockRMSD/examples/data .**
[user@cn3144 ~]$ **cat << EOF > z.py
from pydockrmsd.dockrmsd import PyDockRMSD
import pydockrmsd.hungarian as hungarian
dockrmsd = PyDockRMSD("./data/targets/1a8i/crystal.mol2",
 "./data/targets/1a8i/vina1.mol2")
print(dockrmsd.rmsd)
print(dockrmsd.total\_of\_possible\_mappings)
print(dockrmsd.optimal\_mapping)
print(dockrmsd.error)
print(hungarian("./data/targets/1a8i/crystal.mol2",
 "./data/targets/1a8i/vina1.mol2"))
EOF**
[user@cn3144 ~]$ **python z.py**
0.692337048145839
108.0
Optimal mapping (First file -> Second file, * indicates correspondence is not one-to-one):
C  1 -> C  1
C  2 -> C  2
O  3 -> O 15 *
C  4 -> C  3 *
O  5 -> O 17 *
C  6 -> C  4 *
O  7 -> O 22 *
C  8 -> C  5 *
C  9 -> C 19 *
O 10 -> O 20 *
O 11 -> O  6 *
N 12 -> N  7 *
C 13 -> C  8 *
O 14 -> O  9 *
N 15 -> N 10 *
C 16 -> C 11 *
O 17 -> O 12 *

0.6923022121452281

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. pydockrmsd.sh). For example:



```

#!/bin/bash
set -e
module load pydockrmsd
python my_pydockrmsd_script.py

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch pydockrmsd.sh
```







