

document.querySelector('title').textContent = 'phase on Biowulf';
phase on Biowulf


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



`PHASE` reconstructs haplotypes from population genotype data using
a Bayesian statistical model that considers the decay of LD with distance
due to recombination. Inputs can include biallelic SNPs as well as multi-allelic
loci like SNPs with more than two alleles, HLA allels, or microsatellites.



### References:


* M. Stephens, M. Smith, and P. Donnelly.  *A new statistical method for 
 haplotype reconstruction from population data*. American Journal of Human Genetics
 2001, 68:978-989.
 [Pubmed](http://www.ncbi.nlm.nih.gov/pubmed/11254454) | 
 [PMC](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1275651/) | 
 [Journal](http://www.sciencedirect.com/science/article/pii/S0002929707614244)
* M. Stephens, P. Scheet. *Accounting for Decay of Linkage Disequilibrium in 
 Haplotype Inference and Missing-Data Imputation.* American Journal of Human Genetics
 2005, 76:449-462.
 [Pubmed](http://www.ncbi.nlm.nih.gov/pubmed/15700229) | 
 [PMC](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1196397/) | 
 [Journal](http://www.sciencedirect.com/science/article/pii/S0002929707633412)


Documentation
* [Manual](http://stephenslab.uchicago.edu/assets/software/phase/instruct2.1.pdf)


Important Notes
* Module Name: phase (see [the modules page](/apps/modules.html) for more information)
* Example files in `$PHASE_TEST_DATA`



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int), load
the phase module and analyze a simple test data set running 1000 iterations



```

[user@biowulf]$ **sinteractive --gres=lscratch:5**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3114]$ **module load phase/2.1.1**
[user@cn3114]$ **cd /lscratch/$SLURM\_JOB\_ID**
[user@cn3114]$ **cp $PHASE\_TEST\_DATA/test.inp .**
[user@cn3114]$ **PHASE test.inp test.out 1000**
Reading in data
Reading Positions of loci
Reading individual      3
Finished reading
Computing matrix Q, please wait
Done computing Q
3
5
MSSSM
0 #1
12 0 0 1 3
11 1 1 0 3
0 #2
12 0 1 1 3
12 1 0 0 2
0 #3
12 1 0 1 2
12 0 1 0 13
Resolving with method R
Making List of all possible haplotypes
Method = R
Performing Final Set of Iterations... nearly there!
Performing Burn-in iterations
  50% done
Estimating recom rates
Continuing Burn-in
Performing Main iterations
Writing output to files
Producing Summary, please wait

[user@cn3114]$ **ls -1 test.out\***
test.out
test.out_freqs
test.out_hbg
test.out_monitor
test.out_pairs
test.out_probs
test.out_recom

[user@cn3144]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf]$

```

The input format and options are described in the
[manual](http://stephenslab.uchicago.edu/assets/software/phase/instruct2.1.pdf).



Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. phase.sh), which uses the input file 'phase.in'. For example:



```

#! /bin/bash

module load phase/2.1.1 || exit 1
PHASE input output 1000

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch phase.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. phase.swarm). For example:



```

PHASE -X10 input1 output1
PHASE -X10 input2 output2

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f phase.swarm -g 4 -t 1 -p 2 --module phase/2.1.1
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g #  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t #  Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module phase  Loads the phase module for each subjob in the swarm 
 | |
 | |
 | |








