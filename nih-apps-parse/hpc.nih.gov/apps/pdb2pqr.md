

document.querySelector('title').textContent = "pdb2pqr";
pdb2pqr on Biowulf


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



From the documentation




> 
>  The use of continuum solvation methods such as APBS requires
>  accurate and complete structural data as well as force field parameters
>  such as atomic charges and radii. Unfortunately, the limiting step in
>  continuum electrostatics calculations is often the addition of missing
>  atomic coordinates to molecular structures from the Protein Data Bank and
>  the assignment of parameters to these structures. To address this problem,
>  we have developed PDB2PQR. This software automates many of the common tasks
>  of preparing structures for continuum solvation calculations as well as
>  many other types of biomolecular structure modeling, analysis, and
>  simulation.
> 


### References:


* E. Jurrus, et al. *Improvements to the APBS biomolecular solvation software suite*. 
 Protein Sci. 2018
 [PubMed](https://pubmed.ncbi.nlm.nih.gov/28836357/)  | 
 [PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/pmid/28836357/)  | 
 [Journal](https://onlinelibrary.wiley.com/doi/10.1002/pro.3280)
* T. J. Dolinsky, et al. *PDB2PQR: expanding and upgrading automated preparation of biomolecular structures for molecular simulations*. Nucleic Acids Res 2007
 [PubMed](https://pubmed.ncbi.nlm.nih.gov/17488841/)  | 
 [PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/pmid/17488841/)  | 
 [Journal](https://academic.oup.com/nar/article/35/suppl_2/W522/2920806)


Documentation
* pdb2pqr on [GitHub](https://github.com/Electrostatics/pdb2pqr)
* [Home page](http://www.poissonboltzmann.org)
* [Manual](http://pdb2pqr.readthedocs.io)


Important Notes
* Module Name: pdb2pqr (see [the modules page](/apps/modules.html) for more information)



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --cpus-per-task=2 --mem=2g**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144]$ **pdb2pqr --include-header --ff PARSE --with-ph=7.4 --drop-water 1FAS 1FAS.pqr**
[user@cn3144]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. pdb2pqr.sh), which uses the input file 'pdb2pqr.in'. For example:



```

#!/bin/bash
module load pdb2pqr/3.6.1
pdb2pqr --include-header --ff PARSE --with-ph=7.4 --drop-water 1FAS 1FAS.pqr

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=2 --mem=2g --time 2 --partition=quick pdb2pqr.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. pdb2pqr.swarm). For example:



```

pdb2pqr --include-header --ff PARSE --with-ph=7.4 --drop-water 1FAS 1FAS.pqr
pdb2pqr --include-header --ff PARSE --with-ph=7.4 --drop-water 5EDU 5EDU.pqr
pdb2pqr --include-header --ff PARSE --with-ph=7.4 --drop-water 3C10 3C10.pqr

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f pdb2pqr.swarm -g 2 -t 1 -b3 --time 2 --partition=quick --module pdb2pqr/3.6.1
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g #  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t #  Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module pdb2pqr  Loads the pdb2pqr module for each subjob in the swarm 
 | |
 | |
 | |








