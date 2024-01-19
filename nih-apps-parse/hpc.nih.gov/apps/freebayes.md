

document.querySelector('title').textContent = 'freebayes on Biowulf';
freebayes on Biowulf


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



>  FreeBayes is a Bayesian genetic variant detector designed to find small
> polymorphisms, specifically SNPs (single-nucleotide polymorphisms), indels
> (insertions and deletions), MNPs (multi-nucleotide polymorphisms), and complex
> events (composite insertion and substitution events) smaller than the length of
> a short-read sequencing alignment.
> 
> 
> FreeBayes is haplotype-based, in the sense that it calls variants based on
> the literal sequences of reads aligned to a particular target, not their
> precise alignment. This model is a straightforward generalization of previous
> ones (e.g. PolyBayes, samtools, GATK) which detect or report variants based on
> alignments. This method avoids one of the core problems with alignment-based
> variant detection--- that identical sequences may have multiple possible
> alignments 
> 
> 
> 


### References:


* Erik Garrison, Gabor Marth. *Haplotype-based variant detection from short-read sequencing*.
 2012, [arXiv:1207.3907](https://arxiv.org/abs/1207.3907)


Documentation
* freebayes [GitHub repo](https://github.com/ekg/freebayes)


Important Notes
* Module Name: freebayes (see [the modules page](/apps/modules.html) for more information)
* freebayes can be parallelized with `freebayes-parallel`
* Example data in `$FREEBAYES_TEST_DATA`



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=10g**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load freebayes samtools**
[user@cn3144 ~]$ **freebayes \
 --fasta-reference /fdb/igenomes/Homo\_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa \
 --min-alternate-count 2 \
 --min-alternate-qsum 40 \
 --pvar 0.0001 \
 --use-mapping-quality \
 --site-selection-max-iterations 3 \
 --genotyping-max-iterations 25 \
 /fdb/app\_testdata/bam/hg19/gcat\_set\_053.bam \
 | bgzip -c \
 > test.vcf.gz**

# This takes approximately 72 minutes.

[user@cn3144 ~]$ **bgzip -r test.vcf.gz**
[user@cn3144 ~]$ **ls -lh**
-rw-r--r-- 1 user group 30M Feb 12 14:14 test.vcf.gz
-rw-r--r-- 1 user group 48K Feb 12 14:17 test.vcf.gz.gzi
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. freebayes.sh), which uses the input file 'freebayes.in'. For example:



```

#! /bin/bash

module load freebayes/1.1.0 samtools || exit 1
freebayes --fasta-reference /fdb/igenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa \
    --min-alternate-count 2 \
    --min-alternate-qsum 40 \
    --pvar 0.0001 \
    --use-mapping-quality \
    --site-selection-max-iterations 3 \
    --genotyping-max-iterations 25 \
    secret_sample_1.bam \
  | bgzip -c \
  > secret_sample_1.vcf.gz

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=2 --mem=10g freebayes.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. freebayes.swarm). For example:



```

cd /data/$USER/dir1; freebayes --fasta-reference ../fa/h.sapiens.fasta sample1.bam | bgzip -c > sample1.vcf.gz
cd /data/$USER/dir1; freebayes --fasta-reference ../fa/h.sapiens.fasta sample2.bam | bgzip -c > sample2.vcf.gz
cd /data/$USER/dir1; freebayes --fasta-reference ../fa/h.sapiens.fasta sample3.bam | bgzip -c > sample3.vcf.gz

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f freebayes.swarm -g 10 --module freebayes/1.1.0,samtools
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g #  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t #  Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module freebayes  Loads the freebayes module for each subjob in the swarm 
 | |
 | |
 | |








