

document.querySelector('title').textContent = 'PAINTOR on Biowulf';
PAINTOR on Biowulf


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



PAINTOR is a statistical fine-mapping method that integrates functional genomic data with association strength from potentially multiple populations (or traits) to prioritize variants for follow-up analysis. The software runs on multiple fine-mapping loci and/or populations/traits simultaneously and takes as input the following data for each set of SNPs at a locus
* Summary Association Statistics (Z-scores)
 * Linkage Disequilibrium Matrix/Matrices (Pairwise Pearson correlations coefficients between each SNP)
 * Functional Annotation Matrix (Binary indicator of annotation membership (i.e. if entry {i,k} = 1, then SNP i is a member of annotation K).





### References:


* [Integrating Functional Data to Prioritize Causal Variants in Statistical Fine-Mapping Studies](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1004722). Kichaev et al, PLOS Genetics, 2014
* [Leveraging Functional-Annotation Data in Trans-ethnic Fine-Mapping Studies](https://www.ncbi.nlm.nih.gov/pubmed/?term=26189819). Kichaev and Pasaniuc, Am J Hum Genet, 2015.
* [Improved methods for multi-trait fine mapping of pleiotropic risk loci](https://www.ncbi.nlm.nih.gov/pubmed/?term=27663501). Kichaev et al, Bioinformatics, 2017.


Documentation
* [PAINTOR Main Site](https://github.com/gkichaev/PAINTOR_V3.0)


Important Notes
* Module Name: PAINTOR (see [the modules page](/apps/modules.html) for more information)
 * Singlethreaded
* Example files in /usr/local/apps/PAINTOR/SampleData



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

[user@cn3144 ~]$ **module load PAINTOR**

[user@cn3144 ~]$ **cp -pr /usr/local/apps/PAINTOR/SampleData /data/$USER** 

[user@cn3144 ~]$ **cd /data/$USER**

[user@cn3144 ~]$ **PAINTOR -input SampleData/input.files -in SampleData/ \
 -out SampleData/ -Zhead Zscore -LDname ld -enumerate 2 -annotations DHS**
Running PAINTOR with full enumeration
Maximum number of causals per locus: 2
Proportion of LD variance kept when performing truncated SVD for estimating N*h2g: 0.95
Model annotations: DHS
**********
Reading in files for: Locus1
[....]
Enrichment estimates at iteration 11 is :
 4.15975
-1.88218
Average log bayes factor: 658.648

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. PAINTOR.sh). For example:



```

#!/bin/bash
set -e
module load PAINTOR
cp -pr /usr/local/apps/PAINTOR/SampleData  /data/$USER
PAINTOR -input SampleData/input.files    -in SampleData/  \
      -out SampleData/    -Zhead Zscore -LDname ld    -enumerate 2   -annotations DHS

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--mem=#] PAINTOR.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. PAINTOR.swarm). For example:



```

PAINTOR -input dir1/input.files    -in dir1/  -out dir1  -Zhead Zscore -LDname ld  -enumerate 2 -annotations DHS
PAINTOR -input dir2/input.files    -in dir2/  -out dir1  -Zhead Zscore -LDname ld  -enumerate 2 -annotations DHS
PAINTOR -input dir3/input.files    -in dir3/  -out dir1  -Zhead Zscore -LDname ld  -enumerate 2 -annotations DHS

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f PAINTOR.swarm [-g #] --module PAINTOR
```

where


|  |  |  |  |
| --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | --module PAINTOR Loads the PAINTOR module for each subjob in the swarm 
 | |
 | |








