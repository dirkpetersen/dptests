

document.querySelector('title').textContent = 'RFmix on Biowulf';
RFmix on Biowulf


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



RFMIX is a program to identify the ancestry of genomic segments using random forest discriminative machine learning methods combined with a conditional random field model of the linear chromosome.



### Web site


* [Home page](https://github.com/slowkoni/rfmix)


Documentation
* [RFmix Documentation](https://github.com/slowkoni/rfmix/blob/master/MANUAL.md)


Important Notes
* Module Name: rfmix (see [the modules page](/apps/modules.html) for more information)



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive**
salloc.exe: Pending job allocation 12345678
salloc.exe: job 12345678 queued and waiting for resources
salloc.exe: job 12345678 has been allocated resources
salloc.exe: Granted job allocation 12345678
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn1234 are ready for job

[user@cn1234 ~]$ **module load rfmix**

[user@cn1234 ~]$ **rfmix --help**
RFMIX v2.03-r0 - Local Ancestry and Admixture Inference
(c) 2016, 2017 Mark Koni Hamilton Wright
Bustamante Lab - Stanford University School of Medicine
Based on concepts developed in RFMIX v1 by Brian Keith Maples, et al.

This version is licensed for non-commercial academic research use only
For commercial licensing, please contact cdbadmin@stanford.edu

--- For use in scientific publications please cite original publication ---
Brian Maples, Simon Gravel, Eimear E. Kenny, and Carlos D. Bustamante (2013).
RFMix: A Discriminative Modeling Approach for Rapid and Robust Local-Ancestry
Inference. Am. J. Hum. Genet. 93, 278-288

Summary of command line options - see manual for details

   -f , --query-file= (required)
 VCF file with samples to analyze (required)
 -r , --reference-file= (required)
 VCF file with reference individuals (required)
 -m , --sample-map= (required)
 Reference panel sample population classification map (required)
 -g , --genetic-map= (required)
 Genetic map file (required)
 -o , --output-basename= (required)
 Basename (prefix) for output files (required)
 --chromosome= (required)
 Execute only on specified chromosome (required)

 -c , --crf-spacing=
 Conditional Random Field spacing (# of SNPs)
 -s , --rf-window-size=
 Random forest window size (class estimation window size)
 -w , --crf-weight=
 Weight of observation term relative to transition term in conditional random field
 -G , --generations=
 Average number of generations since expected admixture
 -e , --em-iterations=
 Maximum number of EM iterations
 --reanalyze-reference
 In EM, analyze local ancestry of the reference panel and reclassify it

 -n , --node-size=
 Terminal node size for random forest trees
 -t , --trees=
 Number of tree in random forest to estimate population class probability
 --max-missing=
 Maximum proportion of missing data allowed to include a SNP
 -b , --bootstrap-mode=
 Specify random forest bootstrap mode as integer code (see manual)
 --rf-minimum-snps=
 With genetic sized rf windows, include at least this many SNPs regardless of span
 --analyze-range=
 Physical position range, specified as -, in Mbp (decimal allowed)

 --debug=
 Turn on any debugging output
 --n-threads=
 Force number of simultaneous thread for parallel execution
 --random-seed=
 Seed value for random number generation (integer)
 (maybe specified in hexadecimal by preceeding with 0x), or the string
 "clock" to seed with the current system time.

[user@cn1234 ~]$  **vartrix --bam test\_dna.bam --cell-barcodes dna\_barcodes.tsv --fasta test\_dna.fa --out-matrix test\_dna.out --vcf test\_dna.vcf**

[user@cn1234 ~]$ **exit**
salloc.exe: Relinquishing job allocation 12345678
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. vartrix\_job.sh). For example:



```

#!/bin/bash
set -e
module load vartrix
vartrix --bam test_dna.bam \
   --cell-barcodes dna_barcodes.tsv \
   --fasta test_dna.fa \
   --out-matrix test_dna.out \
   --vcf test_dna.vcf

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] vartrix_job.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. vartrix\_jobs.swarm). For example:



```

vartrix --bam test_dna1.bam \ 
   --cell-barcodes dna_barcodes1.tsv \
   --fasta test_dna1.fa \
   --out-matrix test_dna1.out \
   --vcf test_dna1.vcf
vartrix --bam test_dna2.bam \
   --cell-barcodes dna_barcodes2.tsv \
   --fasta test_dna2.fa \
   --out-matrix test_dna2.out \
   --vcf test_dna2.vcf
vartrix --bam test_dna3.bam \
   --cell-barcodes dna_barcodes3.tsv \
   --fasta test_dna3.fa \
   --out-matrix test_dna3.out \
   --vcf test_dna3.vcf

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f vartrix_jobs.swarm [-g #] [-t #] --module vartrix
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module vartrix Loads the vartrix module for each subjob in the swarm 
 | |
 | |
 | |








