

document.querySelector('title').textContent = 'HipSTR on Biowulf';

 .hl { background-color: #ffff99; }
 code {padding: 1px; background-color: #eeeeee; border: 1px solid #bbbbbbc5;}
 dt {font-weight: bold; margin-top: 5px;}
 dd {padding-left: 2px; border-left: 1px solid #bbbbbb;}

HipSTR on Biowulf


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



 Short tandem repeats (STRs) are highly repetitive genomic sequences comprised of repeated copies of an underlying motif. Prevalent in most organisms' genomes, STRs are of particular interest because they mutate much more rapidly than most other genomic elements. As a result, they're extremely informative for genomic identification, ancestry inference and genealogy.

 Despite their utility, STRs are particularly difficult to genotype. The repetitive sequence responsible for their high mutability also results in frequent alignment errors that can complicate and bias downstream analyses. In addition, PCR stutter errors often result in reads that contain additional or fewer repeat copies than the true underlying genotype.



HipSTR was specifically developed to deal with these errors in the hopes of obtaining more robust STR genotypes. In particular, it accomplishes this by:


* Learning locus-specific PCR stutter models using an EM algorithm 
 * Mining candidate STR alleles from population-scale sequencing data
 * Employing a specialized hidden Markov model to align reads to candidate alleles while accounting for STR artifacts
 * Utilizing phased SNP haplotypes to genotype and phase STRs


### References:


<https://www.nature.com/articles/nmeth.4267>


Documentation
* <https://github.com/tfwillems/HipSTR>


Important Notes
* Module Name: hipstr (see [the modules page](/apps/modules.html) for more information)
* HipSTR is easy to run: `HipSTR -h`
* Example files can be copied from: `/usr/local/apps/hipstr/HipSTR-tutorial`



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive --mem=4g**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144]$ **module load hipstr**
[user@cn3144]$ **cp -r /usr/local/apps/hipstr/HipSTR-tutorial /data/$USER; cd /data/$USER/HipSTR-tutorial**
[user@cn3144]$ **export bam=/data/$USER/HipSTR-tutorial/bams**
[user@cn3144]$ **export fasta=/data/$USER/HipSTR-tutorial/fasta**
[user@cn3144]$ **export tut=/data/$USER/HipSTR-tutorial**
[user@cn3144]$ **HipSTR --bams $bam/ERR194147.bam,$bam/ERR194160.bam,\
 $bam/ERR194161.bam,$bam/SRR826427.bam,\
 $bam/SRR826428.bam,$bam/SRR826448.bam,\
 $bam/SRR826463.bam,$bam/SRR826465.bam,\
 $bam/SRR826467.bam,$bam/SRR826469.bam,\
 $bam/SRR826471.bam,$bam/SRR826473.bam \
 --fasta $fasta/all\_chroms.fa \
 --regions $tut/regions.bed \
 --str-vcf $tut/trio.marshfield.vcf.gz \
 --log $tut/trio.marshfield.log \
 --viz-out $tut/trio.marshfield.viz.gz \
 --min-reads 25 --def-stutter-model**

[user@cn3144]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. hipstr.sh). For example:



```

#!/bin/bash
#SBATCH --job-name=hipstr
#SBATCH --output=hipstr.out
#SBATCH --ntasks=1
#SBATCH --mem=4g
#SBATCH --time=30:00
#SBATCH --partition=quick

set -e
module load hipstr
cp -r /usr/local/apps/hipstr/HipSTR-tutorial /data/$USER
cd /data/$USER/HipSTR-tutorial
HipSTR --bams $bam/ERR194147.bam,$bam/ERR194160.bam,\
              $bam/ERR194161.bam,$bam/SRR826427.bam,\
              $bam/SRR826428.bam,$bam/SRR826448.bam,\
              $bam/SRR826463.bam,$bam/SRR826465.bam,\
              $bam/SRR826467.bam,$bam/SRR826469.bam,\
              $bam/SRR826471.bam,$bam/SRR826473.bam \
        --fasta $fasta/all_chroms.fa \
        --regions $tut/regions.bed \
        --str-vcf $tut/trio.marshfield.vcf.gz \
        --log $tut/trio.marshfield.log \
        --viz-out $tut/trio.marshfield.viz.gz \
        --min-reads 25 --def-stutter-model


```

 Submit the job:

```
sbatch hipstr.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. hipstr.swarm). For example:



```

cd dir1;HipSTR ....... --def-stutter-model
cd dir2;HipSTR ....... --def-stutter-model
cd dir3;HipSTR ....... --def-stutter-model

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f hipstr.swarm -g 8 --module hipstr
```

where


|  |  |
| --- | --- |
| -g #  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | |










