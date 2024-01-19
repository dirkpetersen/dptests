

document.querySelector('title').textContent = 'WHATSHAP on Biowulf';
WHATSHAP on Biowulf


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



WhatsHap is an application for phasing genomic variants using DNA sequencing reads.



### References:


* [WhatsHap: fast and accurate read-based phasing](https://doi.org/10.1101/085050). Marcel Martin, Murray Patterson, Shilpa Garg, Sarah O. Fischer, Nadia Pisanti, Gunnar W. Klau, Alexander Schoenhuth, Tobias Marschall. bioRxiv 085050
* [WhatsHap: Weighted Haplotype Assembly for Future-Generation Sequencing Reads](http://dx.doi.org/10.1089/cmb.2014.0157). Murray Patterson, Tobias Marschall, Nadia Pisanti, Leo van Iersel, Leen Stougie, Gunnar W. Klau, Alexander Schönhuth. Journal of Computational Biology, 22(6), pp. 498-509, 2015.


Documentation
* [Whatshap Main Site](https://whatshap.readthedocs.io)


Important Notes
* Module Name: whatshap (see [the modules page](/apps/modules.html) for more information)



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

[user@cn3144 ~]$ **module load whatshap**

[user@cn3144 ~]$ **cd /data/$USER/WHATSHAP\_TEST**

[user@cn3144 ~]$ **whatshap phase -o phased.vcf input.vcf pacbio\_reads\_30x.bam**
Working on 1 samples from 1 family
======== Working on chromosome 'contig1'
---- Processing individual SAMPLE
Using maximum coverage per sample of 15X
Number of variants skipped due to missing genotypes: 0
Number of remaining heterozygous variants: 177
Reading alignments for sample 'SAMPLE'and detecting alleles ...
WARNING: Sample 'SAMPLE' not found in any BAM/CRAM file.
Found 0 reads covering 0 variants
Kept 0 reads that cover at least two variants each
Reducing coverage to at most 15X by selecting most informative reads ...
Selected 0 reads covering 0 variants
Best-case phasing would result in 0 non-singleton phased blocks (0 in total)
... after read selection: 0 non-singleton phased blocks (0 in total)
Variants covered by at least one phase-informative read in at least one individual after read selection: 0
Phasing 1 sample by solving the MEC problem ...
MEC cost: 0
No. of phased blocks: 0
======== Writing VCF
Done writing VCF
======== Working on chromosome 'contig2'
---- Processing individual SAMPLE
Using maximum coverage per sample of 15X
Number of variants skipped due to missing genotypes: 0
Number of remaining heterozygous variants: 134
Reading alignments for sample 'SAMPLE'and detecting alleles ...
WARNING: Sample 'SAMPLE' not found in any BAM/CRAM file.
Found 0 reads covering 0 variants
Kept 0 reads that cover at least two variants each
Reducing coverage to at most 15X by selecting most informative reads ...
Selected 0 reads covering 0 variants
Best-case phasing would result in 0 non-singleton phased blocks (0 in total)
... after read selection: 0 non-singleton phased blocks (0 in total)
Variants covered by at least one phase-informative read in at least one individual after read selection: 0
Phasing 1 sample by solving the MEC problem ...
MEC cost: 0
No. of phased blocks: 0
======== Writing VCF
Done writing VCF
======== Working on chromosome 'contig3'
---- Processing individual SAMPLE
Using maximum coverage per sample of 15X
Number of variants skipped due to missing genotypes: 0
Number of remaining heterozygous variants: 165
Reading alignments for sample 'SAMPLE'and detecting alleles ...
WARNING: Sample 'SAMPLE' not found in any BAM/CRAM file.
Found 0 reads covering 0 variants
Kept 0 reads that cover at least two variants each
Reducing coverage to at most 15X by selecting most informative reads ...
Selected 0 reads covering 0 variants
Best-case phasing would result in 0 non-singleton phased blocks (0 in total)
... after read selection: 0 non-singleton phased blocks (0 in total)
Variants covered by at least one phase-informative read in at least one individual after read selection: 0
Phasing 1 sample by solving the MEC problem ...
MEC cost: 0
No. of phased blocks: 0
======== Writing VCF
Done writing VCF

== SUMMARY ==
Maximum memory usage: 0.050 GB
Time spent reading BAM/CRAM:                    0.0 s
Time spent parsing VCF:                         0.1 s
Time spent selecting reads:                     0.0 s
Time spent phasing:                             0.0 s
Time spent writing VCF:                         0.1 s
Time spent finding components:                  0.0 s
Time spent on rest:                             0.0 s
Total elapsed time:                             0.2 s

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. whatshap.sh). For example:



```

#!/bin/bash
set -e
module load whatshap
whatshap phase -o /data/$USER/WHATSHAP_TEST/phased.vcf \
                  /data/$USER/WHATSHAP_TEST/input.vcf \
                  /data/$USER/WHATSHAP_TEST/pacbio_reads_30x.bam

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] whatshap.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. whatshap.swarm). For example:



```

whatshap phase -o /data/$USER/WHATSHAP_TEST/phased_A.vcf \
                  /data/$USER/WHATSHAP_TEST/input_A.vcf \
                  /data/$USER/WHATSHAP_TEST/pacbio_reads_A.bam
whatshap phase -o /data/$USER/WHATSHAP_TEST/phased_B.vcf \
                  /data/$USER/WHATSHAP_TEST/input_B.vcf \
                  /data/$USER/WHATSHAP_TEST/pacbio_reads_B.bam

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f whatshap.swarm [-g #] [-t #] --module whatshap
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module whatshap Loads the whatshap module for each subjob in the swarm 
 | |
 | |
 | |








