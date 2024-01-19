

document.querySelector('title').textContent = 'RSD on Biowulf';
RSD on Biowulf


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



Reciprocal Smallest Distance (RSD) is a pairwise orthology algorithm that uses global sequence alignment and maximum likelihood evolutionary distance between sequences to accurately detect orthologs between genomes.



### References:


* Wall, D.P., Fraser, H.B. and Hirsh, A.E.
 **[Detecting putative orthologs](https://www.ncbi.nlm.nih.gov/pubmed/?term=15593400)**
*Bioinformatics (2013), 19:1710-1711.*


Documentation
* RSD Main Site: <https://github.com/todddeluca/reciprocal_smallest_distance>


Important Notes
* Module Name: RSD (see [the modules page](/apps/modules.html) for more information)
* Singlethreaded
* Environment variables set 
	+ RSD\_EXAMPLES* Example files in $RSD\_EXAMPLES



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

[user@cn3144 ~]$ **module load RSD**
[user@cn3144 ~]$ **cp -R $RSD\_EXAMPLES .**
[user@cn3144 ~]$ **rsd\_search -q examples/genomes/Mycoplasma\_genitalium.aa/Mycoplasma\_genitalium.aa \
--subject-genome=examples/genomes/Mycobacterium\_leprae.aa/Mycobacterium\_leprae.aa \
-o Mycoplasma\_genitalium.aa\_Mycobacterium\_leprae.aa\_0.8\_1e-5.orthologs.txt**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. RSD.sh). For example:



```

#!/bin/bash
module load RSD
cp -R $RSD_EXAMPLES .
rsd_search -q examples/genomes/Mycoplasma_genitalium.aa/Mycoplasma_genitalium.aa \
--subject-genome=examples/genomes/Mycobacterium_leprae.aa/Mycobacterium_leprae.aa \
-o Mycoplasma_genitalium.aa_Mycobacterium_leprae.aa_0.8_1e-5.orthologs.txt

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] RSD.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. RSD.swarm). For example:



```

rsd_search -q path_to_org1.aa --subject-genome=path_to_org2.aa -o org1_org2_0.8_1e-5.orthologs.txt
rsd_search -q path_to_org1.aa --subject-genome=path_to_org3.aa -o org1_org3_0.8_1e-5.orthologs.txt
rsd_search -q path_to_org1.aa --subject-genome=path_to_org4.aa -o org1_org4_0.8_1e-5.orthologs.txt
rsd_search -q path_to_org1.aa --subject-genome=path_to_org5.aa -o org1_org5_0.8_1e-5.orthologs.txt

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f RSD.swarm [-g #] [-t #] --module RSD
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module RSD Loads the RSD module for each subjob in the swarm 
 | |
 | |
 | |








