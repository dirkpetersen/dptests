

document.querySelector('title').textContent = 'ZDOCK on Biowulf';
ZDOCK on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
 |



zDOCK performs a full rigid-body search of docking orientations between two proteins. This includes performance optimization and a novel pairwise statistical energy potential.



### References:


* Chen R, Li L, Weng Z
 [**ZDOCK: an initial-stage protein-docking algorithm.**](https://pubmed.ncbi.nlm.nih.gov/12784371/)
*Proteins (2003) Jul 1;52(1):80-7.*


Documentation
* [ZDOCK Main Site](https://zlab.umassmed.edu/zdock/)


Important Notes
* Module Name: zdock (see [the modules page](/apps/modules.html) for more information)
 * **zdock/3.0.2**: single-threaded
 * **zdock/3.0.2\_mpi**: MPI
 * Environment variables set 
	+ ZDOCK\_HOME
	+ OMPI\_MCA\_btl



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

[user@cn3144 ~]$ **module load zdock**

[user@cn3144 ~]$ **zgrep ^ATOM /pdb/pdb/cg/pdb1cgi.ent.gz | grep ' E ' > 1cgi\_r.pdb** 
[user@cn3144 ~]$ **zgrep ^ATOM /pdb/pdb/cg/pdb1cgi.ent.gz | grep ' I ' > 1cgi\_l.pdb** 
[user@cn3144 ~]$ **cp $ZDOCK\_HOME/uniCHARMM .** 
[user@cn3144 ~]$ **mark\_sur 1cgi\_r.pdb 1cgi\_r\_m.pdb** 
[user@cn3144 ~]$ **mark\_sur 1cgi\_l.pdb 1cgi\_l\_m.pdb** 
[user@cn3144 ~]$ **zdock -R 1cgi\_r\_m.pdb -L 1cgi\_l\_m.pdb -o zdock.out** 
[user@cn3144 ~]$ **cp $ZDOCK\_HOME/create\_lig .** 
[user@cn3144 ~]$ **create.pl zdock.out** 

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. zdock.sh). For example:



```

#!/bin/bash
set -e
module load zdock/3.0.2_mpi
cp $ZDOCK_HOME/uniCHARMM .
mark_sur receptor.pdb receptor_m.pdb
mark_sur ligand.pdb ligand_m.pdb
srun --mpi=pmix zdock -R receptor_m.pdb -L ligand_m.pdb -o zdock.out
cp $ZDOCK_HOME/create_lig .
create.pl zdock.out

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --ntasks=16 zdock.sh
```







