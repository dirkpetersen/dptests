

document.querySelector('title').textContent = 'Eigensoft on HPC';
Eigensoft on HPC


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


 **\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*  

 Note, the programs eigenstrat and eigenstratQTL of EIGENSOFT version 2.0  

 have been replaced by smarteigenstrat.perl. Please refer to documentation. 
   

 \*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\***


 


 The EIGENSOFT package combines functionality from population genetics 
 methods ([Patterson 
 et al. 2006](http://genepath.med.harvard.edu/%7Ereich/patterson_eigenanalysis_2006.pdf)) and EIGENSTRAT stratification correction me thod ([Price 
 et al. 2006](http://genepath.med.harvard.edu/%7Ereich/Price%20et%20al.pdf)). The EIGENSTRAT method uses principal components analysis 
 to explicitly model ancestry differences between cases and controls along 
 continuous axe s of variation; the resulting correction is specific to a 
 candidate marker’s variation in frequency across ancestral populations, 
 minimizing spurious associations while maximizing power to detect true associations. 
 The EIGENSOFT package has a built-in plotting script and supports multiple 
 file formats and quantitative phenotypes.Eigensoft was developed at  [Harvard Genetics Department and the Broad Institute](http://genepath.med.harvard.edu/~reich/Software.htm). 


### 


Documentation
* <https://github.com/DReichLab/EIG>



Important Notes
* Module Name: eigensoft (see [the 
 modules page](/apps/modules.html) for more information)
* Example files in /usr/local/apps/eigensoft/version/function 
 where function is one of the categories CONVERTF,EIGENSTRAT,POPGEN.





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

[user@cn3144 ~]$ **module load eigensoft**
[user@cn3144 ~]$ **cp -rp /usr/local/apps/eigensoft/6.1.4/CONVERTF /data/$USER**[user@cn3144 ~]$ **cd /data/$USER/CONVERTF**[user@cn3144 ~]$ **perl example.perl**
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
module load eigensoft
eigensoft commands
```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch batch.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. job.swarm). For example:



```

cd dir1; eigensoft commands
cd dir2; egiensoft commands
cd dir3; eigensoft commands

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f job.swarm --module eigensoft
```

where
 

|  |  |
| --- | --- |
| -g *#*  | Number of Gigabytes of memory required for each process 
 (1 line in the swarm command file)  |
| --module  | Loads the module for each subjob in the swarm  |




