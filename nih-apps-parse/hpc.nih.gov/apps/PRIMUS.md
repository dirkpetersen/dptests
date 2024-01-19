

document.querySelector('title').textContent = 'PRIMUS on Biowulf';
PRIMUS on Biowulf


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



PRIMUS stands for Pedigree Reconstruction (PR) and Identification of a Maximum Unrelated Set (IMUS).

> 
>  The IMUS method is an algorithm adapted from graph theory that always identifies the maximum set of unrelated individuals
>  in any dataset, and allows weighting parameters to be utilized in unrelated sample selection. PRIMUS reads in user-generated
>  IBD estimates and outputs the maximum possible set of unrelated individuals, given a specified threshold of relatedness.
>  Additional information for preferential selection of individuals may also be utilized.
> 



> 
>  The PR algorithm is a method to reconstruct pedigrees within a genetic dataset. PRIMUS can verify expected pedigree structures
>  from genetic data, and it can identify and incorporate novel, cryptic relationships into pedigrees.
> 





### References:


* Jeffrey Staples, Dandi Qiao, Michael H. Cho, Edwin K. Silverman, University of Washington Center for Mendelian Genomics,
 Deborah A. Nickerson, and Jennifer E. Below. 
 [**PRIMUS: Rapid reconstruction of pedigrees from genome-wide estimates of identity by descent.**](https://doi.org/10.1016/j.ajhg.2014.10.005)
*The American Journal of Human Genetics, Volume 95, Issue 5, 2014, Pages 553-564*


Documentation
* [PRIMUS Main Site](https://primus.gs.washington.edu/primusweb/index.html)


Important Notes
* Module Name: PRIMUS (see [the modules page](/apps/modules.html) for more information)
* Singlethreaded
* For example data, see the environment variable, PRIMUS\_TEST\_DATA



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

[user@cn3144 ~]$ **module load PRIMUS**
[+] Loading plink  2.3-alpha
[+] Loading PRIMUS  1.9.0
[user@cn3144 ~]$ **cp $PRIMUS\_TEST\_DATA/complete.genome .**
[user@cn3144 ~]$ **run\_PRIMUS.pl --plink complete.genome**

FILES AND COLUMNS
LOG FILE: complete.genome_PRIMUS/PRIMUS_output.log
IBD file: complete.genome (FID1=1; IID1=2; FID2=3; IID2=4; IBD0=7; IBD1=8; IBD2=9; PI_HAT/RELATEDNESS=10)
Dataset results dir: complete.genome_PRIMUS
Age file: none
Sex file: none
Affection file: none
Trait weighting:
        size (size)

SETTINGS
Get PLINK IBD ESTIMATES with prePRIMUS: 0
Automatic reference population selection: 1
Verbose: 1
Relatedness threshold: 0.09375
Initial likelihood cutoff: 0.3
Max generations: none
Max generational mating gap: 0
Get max unrelated set: 1
Reconstruct pedigrees: 1
Relatedness_file: complete.genome
Threshold: 0.09375
Selection criteria are based on the following:
        size (size)

IDENTFYING FAMILY NETWORKS IN DATA
Writing network files to complete.genome_PRIMUS/
Loading data...
done.
done.

IDENTIFYING A MAXIMUM UNRELATED SET
Checking for large networks...
done.
# of family networks: 1
Writing out unrelated set
done.
Testing alternative methods...
done.
unrelated_file: complete.genome_maximum_independent_set
unrelated_set size: 6
RECONSTRUCTING complete.genome_network1
Output directory: complete.genome_PRIMUS/complete.genome_network1
Use mito non-match: 0
Use mito match: 0
Use Y non-match: 0
Use Y match: 0
Entering Resolve PC trios. # of possible pedigrees: 1
Entering Phase 1. # of possible pedigrees: 1
Entering Phase 2. # of possible pedigrees: 1
Entering Phase 3. # of possible pedigrees: 1
networks pre-prune: 1
networks post-prune: 1
Writing summary file
Writing .fam file for complete.genome_network1_1
Writing dataset Summary file complete.genome_PRIMUS/Summary_complete.genome.txt
done.
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. PRIMUS.sh). For example:



```

#!/bin/bash
set -e
module load PRIMUS
run_PRIMUS.pl --plink input.genome

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--mem=#] PRIMUS.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. PRIMUS.swarm). For example:



```

run_PRIMUS.pl --plink input1.genome
run_PRIMUS.pl --plink input2.genome
run_PRIMUS.pl --plink input3.genome
run_PRIMUS.pl --plink input4.genome

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f PRIMUS.swarm [-g #] --module PRIMUS
```

where


|  |  |  |  |
| --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | --module PRIMUS Loads the PRIMUS module for each subjob in the swarm 
 | |
 | |








