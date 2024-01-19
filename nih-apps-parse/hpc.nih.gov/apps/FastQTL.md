

document.querySelector('title').textContent = 'Fast and efficient QTL mapping';
**Fast and efficient QTL mapping**


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



In order to discover quantitative trait loci (QTLs), 
multi-dimensional genomic datasets combining
DNA-seq and ChiP-/RNA-seq require methods that rapidly correlate tens of thousands of
molecular phenotypes with millions of genetic variants while appropriately controlling for multiple testing. FastQTL implements a popular cis-QTL mapping strategy
in a user- and cluster-friendly tool. FastQTL also proposes an efficient 
permutation procedure to control for multiple testing. 



### References:


* Halit Ongen, Alfonso Buil, Andrew Anand Brown,
Emmanouil T. Dermitzakis and Olivier Delaneau,
"Fast and efficient QTL mapper for thousands of
molecular phenotypes",  *Bioinformatics*, **32**(10), 2016, 1479–1485


Documentation
* [FastQTL SourceForge page](http://fastqtl.sourceforge.net/)


Important Notes
* Module Name: FastQTL (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Fast; with a permutation scheme relying on Beta approximation.
* Flexible; association testing is implemented w/o qualitative/quantitative covariates.
* User-friendly; only standard file formats are used and easy-to-use options implemented.
* Cluster-friendly; parallelization is easy to set up for large deployment on compute clusters.
* Unusual environment variables set
	+ **FASTQTL\_HOME**  installation directory
	+ **FASTQTL\_BIN**       executable directory
	+ **FASTQTL\_SRC**      source code directory
	+ **FASTQTL\_TEST**    sample data directory
Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=4g**
[user@@cn3316 ~]$ **module load FastQTL**
[user@@cn3316 ~]$ **cp ${FASTQTL\_TEST}/\* .** 
[user@@cn3316 ~]$ **fastQTL -V genotypes.vcf.gz -B phenotypes.bed.gz -O res -L res.log --chunk 1 10** 
Fast QTL
  * Authors : Olivier DELANEAU, Halit ONGEN, Alfonso BUIL & Manolis DERMITZAKIS
  * Contact : olivier.delaneau@gmail.com
  * Webpage : http://fastqtl.sourceforge.net/
  * Version : v2.0

Perform nominal analysis (used to get raw p-values of association)
  * Using p-value threshold = 1.0000000000
  * Random number generator is seeded with 1525965118
  * Considering variants within 1e+06 bp of the MPs
  * Chunk processed 1 / 10

Scanning phenotype data in [phenotypes.bed.gz]
  * 364 phenotypes

Reading phenotype data in [phenotypes.bed.gz]
  * region = 22:17517460-20748405
  * 373 samples included
  * 45 phenotypes included
Reading genotype data in [genotypes.vcf.gz] in VCF format
  * region = 22:16517460-21748405
  * 373 samples included
  * 18602 sites included

Imputing missing genotypes

Imputing missing phenotypes

Initialize covariate 

Processing gene [ENSG00000237438.1]
  * Number of variants in cis = 7966
  * Progress = 2.2%

Processing gene [ENSG00000177663.8]
  * Number of variants in cis = 8165
  * Progress = 4.4%

Processing gene [ENSG00000183307.3]
  * Number of variants in cis = 8279
  * Progress = 6.7%

Processing gene [ENSG00000069998.8]
  * Number of variants in cis = 8466
  * Progress = 8.9%
...
...
Processing gene [ENSG00000206176.5]
  * Number of variants in cis = 7031
  * Progress = 93.3%

Processing gene [ENSG00000196622.5]
  * Number of variants in cis = 7212
  * Progress = 95.6%

Processing gene [ENSG00000161133.12]
  * Number of variants in cis = 6103
  * Progress = 97.8%

Processing gene [ENSG00000185252.13]
  * Number of variants in cis = 6090
  * Progress = 100.0%

Running time: 97 seconds

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. FastQTL.sh). For example:



```

#!/bin/bash
module load FastQTL
fastQTL -V genotypes.vcf.gz -B  phenotypes.bed.gz -O res -L res.log --chunk 7 10

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] FastQTL.sh
```
