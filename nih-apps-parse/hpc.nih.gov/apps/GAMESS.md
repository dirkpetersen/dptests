

document.querySelector('title').textContent = 'GAMESS on Biowulf';
GAMESS on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
 |


The General Atomic and Molecular Electronic Structure System ([GAMESS](http://www.msg.ameslab.gov/GAMESS/)) is a general ab initio quantum
chemistry package. GAMESS is maintained by the members of the Gordon research
group at Iowa State University.


### References:


* M.W.Schmidt, K.K.Baldridge, J.A.Boatz, S.T.Elbert, M.S.Gordon, J.H.Jensen, S.Koseki, N.Matsunaga, K.A.Nguyen, S.J.Su, T.L.Windus, M.Dupuis, J.A.Montgomery
 **General Atomic and Molecular Electronic Structure System**
*J. Comput. Chem. 14, 1347-1363 (1993).*


Documentation
* <https://www.msg.chem.iastate.edu/gamess/documentation.html>


Important Notes
* Module Name: GAMESS (see [the modules page](/apps/modules.html) for more information)
* Multithreaded
* environment variables set 
	+ GAMESS\_HOME* Example files in $GAMESS\_HOME/tests



The arguments to the **rungms** command are:


*rungms inp [version] [nprocs]*


where,


* *inp* - input file name (.inp extension assumed)
* *version* - GAMESS version number (optional, defaults to system
default which is always 00)
* *nprocs* - total number of processes to run (optional, defaults to 1)


Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --cpus-per-task=16 --mem-per-cpu=1024**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ module load GAMESS
[user@cn3144 ~]$ rungms mystuff 00 16

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. GAMESS.sh). For example:



```

#!/bin/bash
module load GAMESS
rungms scf_44 00 $SLURM_CPUS_PER_TASK
```


Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] GAMESS.sh
```







