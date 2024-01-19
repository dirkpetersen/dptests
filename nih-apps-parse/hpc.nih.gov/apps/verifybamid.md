

document.querySelector('title').textContent = 'verifybamid on Biowulf';
verifybamid on Biowulf



|  |
| --- |
| 
Quick Links
[Interactive job on Biowulf](#int)
[Batch job on Biowulf](#serial)
[Swarm of jobs](#swarm)
 |



Description
 verifyBamID verifies whether the reads in particular bam file match
previously known genotypes for an individual (or a group of individuals). It
also checks for sample swaps and cross-contaimination between samples.


There may be multiple versions of verifybamid available. An easy way of selecting the
version is to use [modules](/apps/modules.html). To see the modules
available, type



```

module avail verifybamid 

```

To select a module use



```

module load verifybamid/[version]

```

where `[version]` is the version of choice.



verifybamid is a multithreaded application. Make sure to match the
number of cpus requested with the number of threads.


### Environment variables set


* `$PATH`


### References


* G. Jun, M. Flickinger, K. N. Hetrick, Kurt, J. M. Romm, K. F. Doheny, G. Abecasis, M. Boehnke,and H. M. Kang. *Detecting and Estimating Contamination of Human DNA Samples in Sequencing and Array-Based Genotype Data*. American Journal of Human Genetics 2012, 91:839-848.
 [PubMed](https://www.ncbi.nlm.nih.gov/pubmed/23103226)  | 
 [PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3487130/)  | 
 [Journal](http://www.sciencedirect.com/science/article/pii/S0002929712004788)


### Documentation


* [Manual](http://genome.sph.umich.edu/wiki/VerifyBamID)
* [GitHub](https://github.com/statgen/verifyBamID)




Interactive job on Biowulf
Allocate an interactive session with [sinteractive](/docs/userguide.html#int)
and use as follows



```

biowulf$ **sinteractive** 
node$ **module load verifybamid**
verifyBamID 1.1.3 -- verify identity and purity of sequence data
(c) 2010-2014 Hyun Min Kang, Goo Jun, and Goncalo Abecasis


Available Options
                             Input Files : --vcf [], --bam [], --bai [],
                                           --subset [], --smID []
                    VCF analysis options : --genoError [1.0e-03],
                                           --minAF [0.01],
                                           --minCallRate [0.50]
   Individuals to compare with chip data : --site, --self, --best
          Chip-free optimization options : --free-none, --free-mix [ON],
                                           --free-refBias, --free-full
          With-chip optimization options : --chip-none, --chip-mix [ON],
                                           --chip-refBias, --chip-full
                    BAM analysis options : --ignoreRG, --ignoreOverlapPair,
                                           --noEOF, --precise, --minMapQ [10],
                                           --maxDepth [20], --minQ [13],
                                           --maxQ [40], --grid [0.05]
                 Modeling Reference Bias : --refRef [1.00], --refHet [0.50],
                                           --refAlt [0.00]
                          Output options : --out [], --verbose
                               PhoneHome : --noPhoneHome,
                                           --phoneHomeThinning [50]

node$ **verifyBamID --vcf input1.vcf --bam input1.bam --out output1 --verbose --ignoreRG**
node$ **exit**
biowulf$

```



Batch job on Biowulf
Create a batch script similar to the following example:



```

#! /bin/bash
# this file is verifybamid.sh

module load verifybamid || exit 1
verifyBamID --vcf input1.vcf --bam input1.bam --out output1 \
    --verbose --ignoreRG --noPhoneHome

```

Submit to the queue with [sbatch](/docs/userguide.html):



```

biowulf$ **sbatch verifybamid.sh**

```



Swarm of jobs on Biowulf
Create a swarm command file similar to the following example:



```

# this file is verifybamid.swarm
verifyBamID --vcf input1.vcf --bam input1.bam --out output1 --verbose --ignoreRG --noPhoneHome
verifyBamID --vcf input2.vcf --bam input2.bam --out output2 --verbose --ignoreRG --noPhoneHome
verifyBamID --vcf input3.vcf --bam input3.bam --out output3 --verbose --ignoreRG --noPhoneHome

```

And submit to the queue with [swarm](/apps/swarm.html)



```

biowulf$ **swarm -f verifybamid.swarm**

```



