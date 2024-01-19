

document.querySelector('title').textContent = ' GCTA on Biowulf ';

GCTA on Biowulf 



|  |
| --- |
| 
Quick Links
[Batch job on Biowulf](#serial)
[Swarm of jobs](#swarm)
[Interactive job on Biowulf](#int)
 |



Description
GCTA (Genome-wide Complex Trait Analysis) was originally designed to
estimate the proportion of phenotypic variance explained by genome- or
chromosome-wide SNPs for complex traits (the GREML method), and has
subsequently extended for many other analyses to better understand the genetic
architecture of complex traits. GCTA currently supports the following
functionalities: 


* Estimate the genetic relationship from genome-wide SNPs
* Estimate the inbreeding coefficient from genome-wide SNPs
* Estimate the variance explained by all the autosomal SNPs
* Partition the genetic variance onto individual chromosomes
* Estimate the genetic variance associated with the X-chromosome
* Test the effect of dosage compensation on genetic variance on the X-chromosome
* Predict the genome-wide additive genetic effects for individual subjects and for individual SNPs
* Estimate the LD structure encompassing a list of target SNPs
* Simulate GWAS data based upon the observed genotype data
* Convert Illumina raw genotype data into PLINK format
* Conditional & joint analysis of GWAS summary statistics without individual level genotype data
* Estimating the genetic correlation between two traits (diseases) using SNP data
* Mixed linear model association analysis


There may be multiple versions of GCTA available. An easy way of selecting the
version is to use [modules](/apps/modules.html). To see the modules
available, type



```

module avail GCTA 

```

To select a module use



```

module load GCTA/[version]

```

where `[version]` is the version of choice.



GCTA is a multithreaded application. Make sure to match the
number of cpus requested with the number of threads.


### Environment variables set


* `$PATH`
* `$GCTAHOME`: location of test data


### References


* J. Yang J, S. H. Lee, M. E. Goddard and P. M. Visscher PM. 
 *GCTA: a tool for Genome-wide Complex Trait Analysis*. 
 Am J Hum Genet. 2011, 76-82.
 [Pubmed](https://www.ncbi.nlm.nih.gov/pubmed/21167468)  | 
 [PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3014363/)  | 
 [Journal](http://www.sciencedirect.com/science/article/pii/S0002929710005987)


### Documentation


* [Home page](http://cnsgenomics.com/software/gcta/index.html)
* [Manual](http://cnsgenomics.com/software/gcta/GCTA_UserManual_v1.24.pdf)
* [Forum](http://gcta.freeforums.net/)




Batch job on Biowulf
Create a batch script similar to the following example:



```

#! /bin/bash
# this file is GCTA.batch
module load GCTA
gcta --bfile test --make-grm --out test
```


Submit to the queue with [sbatch](/docs/userguide.html):



```

biowulf$ **sbatch --cpus-per-task=1 --mem-per-cpus=2g GCTA.batch**

```



Swarm of jobs on Biowulf
Create a swarm command file similar to the following example:



```

# this file is GCTA.swarm
gcta --bfile test2 --make-grm --out test2
gcta --bfile test3 --make-grm --out test3
gcta --bfile test4 --make-grm --out test4

```

And submit to the queue with [swarm](/apps/swarm.html)



```

biowulf$ **swarm -f GCTA.swarm --module GCTA -g 2**

```



Interactive job on Biowulf
Allocate an interactive session with [sinteractive](/docs/userguide.html#int)
and use as described [above](#helix)



```

biowulf$ **sinteractive --mem=4g --cpus-per-task=2**
node$ **module load GCTA**
node$ **cp $GCTAHOME/test.\* .**
node$ **gcta --bfile test --make-grm --out test**
node$ **exit**
biowulf$

```





