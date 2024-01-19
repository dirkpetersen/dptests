

document.querySelector('title').textContent = 'Merlin on HPC';
Merlin on HPC  MERLIN uses sparse trees to represent gene flow in pedigrees and is one of 
 the fastest pedigree analysis packages around ([Abecasis et al, 2002](https://www.ncbi.nlm.nih.gov/pubmed/?term=11731797)). 


### 


Documentation
* <http://www.sph.umich.edu/csg/abecasis/Merlin/tour/linkage.html>



Important Notes
* Module Name: merlin (see [the modules 
 page](/apps/modules.html) for more information)
* Example files: /usr/local/apps/merlin/version/examples/





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

[user@cn3144 ~]$ **module load merlin**
[user@cn3144 ~]$ **merlin -d x.dat -p x.ped -m x.map**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```




Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. batch.sh). For example:



```

#!/bin/bash
set -e
module load merlin
merlin -d x.dat -p x.ped -m x.map
```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch batch.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. job.swarm). For example:



```

cd dir1; merlin -d x.dat -p x.ped -m x.map
cd dir2; merlin -d x.dat -p x.ped -m x.map
cd dir3; merlin -d x.dat -p x.ped -m x.map

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f job.swarm --module merlin
```



