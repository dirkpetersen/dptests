

document.querySelector('title').textContent = 'XHMM on Biowulf';
XHMM on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
 |


The XHMM (eXome-Hidden Markov Model) C++ software suite was written to call copy number variation (CNV) from next-generation sequencing projects, where exome capture was used (or targeted sequencing, more generally).


XHMM uses principal component analysis (PCA) normalization and a hidden Markov model (HMM) to detect and genotype copy number variation (CNV) from normalized read-depth data from targeted sequencing experiments.


XHMM was explicitly designed to be used with targeted exome sequencing at high coverage (at least 60x - 100x) using Illumina HiSeq (or similar) sequencing of at least ~50 samples. However, no part of XHMM explicitly requires these particular experimental conditions, just high coverage of genomic regions for many samples.


### References:


* Menachem Fromer, Jennifer L. Moran, Kimberly Chambert, Eric Banks, Sarah E. Bergen, Douglas M. Ruderfer, Robert E. Handsaker, Steven A. McCarroll, Michael C. O'Donovan, Michael J. Owen, George Kirov, Patrick F. Sullivan, Christina M. Hultman, Pamela Sklar, and Shaun M. Purcell.
 [**Discovery and statistical genotyping of copy-number variation from whole-exome sequencing depth.**](http://www.cell.com/AJHG/abstract/S0002-9297%2812%2900417-X)
*American Journal of Human Genetics, 91:597-607, Oct 2012.*
* Christopher S. Poultney, Arthur P. Goldberg, Elodie Drapeau, Yan Kou, Hala Harony-Nicolas, Yuji Kajiwara, Silvia De Rubeis, Simon Durand, Christine Stevens, Karola Rehnstrom, Aarno Palotie, Mark J. Daly, Avi Ma'ayan, Menachem Fromer, and Joseph D. Buxbaum.
 [**Identification of small exonic CNV from whole-exome sequence data and application to autism spectrum disorder.**](http://www.sciencedirect.com/science/article/pii/S000292971300414X)
*American Journal of Human Genetics, 93(4):607-619, 2013.*
* Menachem Fromer and Shaun M. Purcell.
 [**Using XHMM software to detect copy number variation in whole-exome sequencing data.**](http://dx.doi.org/10.1002/0471142905.hg0723s81)
*In Current Protocols in Human Genetics. John Wiley and Sons, Inc., 2014.*


Documentation
* [XHMM Tutorial](http://atgu.mgh.harvard.edu/xhmm/tutorial.shtml)


Important Notes
* Module Name: XHMM (see [the modules page](/apps/modules.html) for more information)


A **params.txt** will need to be created. Here is an example:



```
1e-8	6	70	-3	1.00	0	1.00	3	1.00
```

A parameters file consists of the following 9 values (tab-delimited):


1. Exome-wide CNV rate
2. Mean number of targets in CNV
3. Mean distance between targets within CNV (in KB)
4. Mean of DELETION z-score distribution
5. Standard deviation of DELETION z-score distribution
6. Mean of DIPLOID z-score distribution
7. Standard deviation of DIPLOID z-score distribution
8. Mean of DUPLICATION z-score distribution
9. Standard deviation of DUPLICATION z-score distribution


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

[user@cn3144 ~]$ **module load XHMM**
[user@cn3144 ~]$ **xhmm -p params.txt**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. XHMM.sh). For example:



```

#!/bin/bash
module load XHMM
xhmm --mergeGATKdepths -o ./DATA.RD.txt \
--GATKdepths group1.DATA.sample_interval_summary \
--GATKdepths group2.DATA.sample_interval_summary \
--GATKdepths group3.DATA.sample_interval_summary

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] XHMM.sh
```





