

document.querySelector('title').textContent = 'LONGSHOT on Biowulf';
LONGSHOT on Biowulf


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



Longshot is a variant calling tool for diploid genomes using long error prone reads such as Pacific Biosciences (PacBio) SMRT and Oxford Nanopore Technologies (ONT). It takes as input an aligned BAM file and outputs a phased VCF file with variants and haplotype information. It can also output haplotype-separated BAM files that can be used for downstream analysis. Currently, it only calls single nucleotide variants (SNVs).



### References:


* [Longshot: accurate variant calling in diploid genomes using single-molecule long read sequencing](https://www.biorxiv.org/content/10.1101/564443v1). Edge, P. and Bansal, V., 2019. bioRxiv, p.564443.


Documentation
* [Longshot Main Site](https://github.com/pjedge/longshot)


Important Notes
* Module Name: longshot (see [the modules page](/apps/modules.html) for more information)



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

[user@cn3144 ~]$ **module load longshot**

[user@cn3144 ~]$ **cd /data/$USER/LONGSHOT\_TEST**

[user@cn3144 ~]$ **ls**
genome.fa  genome.fa.fai  pacbio_reads_30x.bam	pacbio_reads_30x.bam.bai

[user@cn3144 ~]$ **longshot --bam pacbio\_reads\_30x.bam --ref genome.fa --out longshot\_output.vcf**

2019-11-26 09:58:25 Min read coverage set to 6.
2019-11-26 09:58:25 Max read coverage set to 8000.
2019-11-26 09:58:25 Estimating alignment parameters...
2019-11-26 09:58:26 Done estimating alignment parameters.

                    Transition Probabilities:
                    match -> match:          0.871
                    match -> insertion:      0.099
                    match -> deletion:       0.030
                    deletion -> match:       0.964
                    deletion -> deletion:    0.036
                    insertion -> match:      0.894
                    insertion -> insertion:  0.106

                    Emission Probabilities:
                    match (equal):           0.978
                    match (not equal):       0.007
                    insertion:               1.000
                    deletion:                1.000
                    GENOTYPE PRIORS:
                    REF G1/G2 PROB
                    C D/I 0
                    G A/A 0.00016666692910805806
                    G D/I 0
                    T T/T 0.9985001541646916
                    A C/D 0
                    A A/T 0.00033333385823389856
[...]
2019-11-26 09:58:26 Calling potential SNVs using pileup...
2019-11-26 09:58:27 748 potential SNVs identified.
2019-11-26 09:58:27 Generating haplotype fragments from reads...
2019-11-26 09:58:27    10% of variants processed...
2019-11-26 09:58:28    20% of variants processed...
2019-11-26 09:58:28    30% of variants processed...
2019-11-26 09:58:28    40% of variants processed...
2019-11-26 09:58:28    50% of variants processed...
2019-11-26 09:58:28    60% of variants processed...
2019-11-26 09:58:28    70% of variants processed...
2019-11-26 09:58:28    80% of variants processed...
2019-11-26 09:58:28    90% of variants processed...
2019-11-26 09:58:29    100% of variants processed.
2019-11-26 09:58:29 Calling initial genotypes using pair-HMM realignment...
2019-11-26 09:58:29 Iteratively assembling haplotypes and refining genotypes...
2019-11-26 09:58:29    Round 1 of haplotype assembly...
2019-11-26 09:58:29    (Before HapCUT2) Total phased heterozygous SNVs: 468  Total likelihood (phred): 211782.83
2019-11-26 09:58:31    (After HapCUT2)  Total phased heterozygous SNVs: 468  Total likelihood (phred): 39775.84
2019-11-26 09:58:31    (After Greedy)   Total phased heterozygous SNVs: 468  Total likelihood (phred): 39507.12
2019-11-26 09:58:31    Round 2 of haplotype assembly...
2019-11-26 09:58:31    (Before HapCUT2) Total phased heterozygous SNVs: 476  Total likelihood (phred): 39507.12
2019-11-26 09:58:32    (After HapCUT2)  Total phased heterozygous SNVs: 476  Total likelihood (phred): 39507.12
2019-11-26 09:58:33    (After Greedy)   Total phased heterozygous SNVs: 476  Total likelihood (phred): 39507.12
2019-11-26 09:58:33 Printing VCF file...

[user@cn3144 ~]$ **ls**
genome.fa  genome.fa.fai  longshot_output.vcf  pacbio_reads_30x.bam  pacbio_reads_30x.bam.bai

[user@cn3144 ~]$ **head longshot\_output.vcf**
##fileformat=VCFv4.2
##source=Longshot v0.3.5
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth of reads passing MAPQ filter">
##INFO=<ID=AC,Number=R,Type=Integer,Description="Number of Observations of Each Allele">
##INFO=<ID=AM,Number=1,Type=Integer,Description="Number of Ambiguous Allele Observations">
##INFO=<ID=MC,Number=1,Type=Integer,Description="Minimum Error Correction (MEC) for this single variant">
##INFO=<ID=MF,Number=1,Type=Float,Description="Minimum Error Correction (MEC) Fraction for this variant.">
##INFO=<ID=MB,Number=1,Type=Float,Description="Minimum Error Correction (MEC) Fraction for this variant's haplotype block.">
##INFO=<ID=AQ,Number=1,Type=Float,Description="Mean Allele Quality value (PHRED-scaled).">
##INFO=<ID=GM,Number=1,Type=Integer,Description="Phased genotype matches unphased genotype (boolean).">

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226

[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. longshot.sh). For example:



```

#!/bin/bash
set -e
module load longshot
longshot --bam /data/$USER/LONGSHOT_TEST/pacbio_reads_30x.bam \
         --ref /data/$USER/LONGSHOT_TEST/genome.fa \
         --out /data/$USER/LONGSHOT_TEST/longshot_output.vcf

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] longshot.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. longshot.swarm). For example:



```

longshot --bam /data/$USER/LONGSHOT_TEST/pacbio_reads_1.bam \
         --ref /data/$USER/LONGSHOT_TEST/genome.fa \
         --out /data/$USER/LONGSHOT_TEST/output_1.vcf
longshot -A -r chr1 \
         --bam /data/$USER/LONGSHOT_TEST/pacbio_reads_2.bam \
         --ref /data/$USER/LONGSHOT_TEST/genome.fa \
         --out /data/$USER/LONGSHOT_TEST/output_2.vcf

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f longshot.swarm [-g #] [-t #] --module longshot
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module longshot Loads the longshot module for each subjob in the swarm 
 | |
 | |
 | |








