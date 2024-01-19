

document.querySelector('title').textContent = 'GLNEXUS on Biowulf';
GLNEXUS on Biowulf


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



Scalable gVCF merging and joint variant calling for population sequencing projects.



### References:


* [GLnexus: joint variant calling for large cohort sequencing.](https://doi.org/10.1101/343970) 
Michael F. Lin, Ohad Rodeh, John Penn, Xiaodong Bai, Jeffrey G. Reid, Olga Krasheninina, William J. Salerno.
bioRxiv 343970.


Documentation
* [GLnexus Main Site](https://github.com/dnanexus-rnd/GLnexus)


Important Notes
* Module Name: glnexus (see [the modules page](/apps/modules.html) for more information)
 * Note that as of version 1.4.1, the glnexus command has been changed to glnexus\_cli



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
[user@cn3144 ~]$ **cp ${GLNEXUS\_TEST\_DATA}/\* .**
[user@cn3144 ~]$ **tar -xvf dv\_platinum6\_chr21\_gvcf.tar**
dv_platinum6_chr21_gvcf/
dv_platinum6_chr21_gvcf/NA12890.chr21.gvcf.gz
dv_platinum6_chr21_gvcf/NA12892.chr21.gvcf.gz
dv_platinum6_chr21_gvcf/NA12891.chr21.gvcf.gz
dv_platinum6_chr21_gvcf/NA12889.chr21.gvcf.gz
dv_platinum6_chr21_gvcf/NA12877.chr21.gvcf.gz
dv_platinum6_chr21_gvcf/NA12878.chr21.gvcf.gz

[user@cn3144 ~]$ **rm dv\_platinum6\_chr21\_gvcf.tar**
[user@cn3144 ~]$ **ls**
dv_platinum6_chr21_gvcf

[user@cn3144 ~]$ **module load glnexus**
[user@cn3144 ~]$ **echo -e "chr21\t0\t48129895" > hg19\_chr21.bed**
[user@cn3144 ~]$ **glnexus\_cli --config DeepVariant --bed hg19\_chr21.bed dv\_platinum6\_chr21\_gvcf/\*.gvcf.gz > dv\_platinum6\_chr21.bcf**
[1485789] [2023-09-08 12:31:07.444] [GLnexus] [info] glnexus_cli release v1.4.1-0-g68e25e5 Aug 13 2021
[1485789] [2023-09-08 12:31:07.445] [GLnexus] [info] detected jemalloc 5.2.1-0-gea6b3e973b477b8061e0076bb257dbd7f3faa756
[1485789] [2023-09-08 12:31:07.446] [GLnexus] [info] Loading config preset DeepVariant
[1485789] [2023-09-08 12:31:07.449] [GLnexus] [info] config:
unifier_config:
  drop_filtered: false
  min_allele_copy_number: 1
  min_AQ1: 10
  min_AQ2: 10
  min_GQ: 0
  max_alleles_per_site: 32
  monoallelic_sites_for_lost_alleles: true
  preference: common
genotyper_config:
  revise_genotypes: true
  min_assumed_allele_frequency: 9.99999975e-05
  snv_prior_calibration: 0.600000024
  indel_prior_calibration: 0.449999988
  required_dp: 0
  allow_partial_data: true
  allele_dp_format: AD
  ref_dp_format: MIN_DP
  output_residuals: false
  more_PL: true
  squeeze: false
  trim_uncalled_alleles: true
  top_two_half_calls: false
  output_format: BCF
  liftover_fields:

[...]

[1485789] [2023-09-08 12:31:07.605] [GLnexus] [info] db_get_contigs GLnexus.DB
[1485789] [2023-09-08 12:31:07.674] [GLnexus] [info] Beginning bulk load with no range filter.
[1485789] [2023-09-08 12:31:10.919] [GLnexus] [info] Loaded 6 datasets with 6 samples; 239726128 bytes in 2572592 BCF records (10 duplicate) in 7062 buckets. Bucket max 551480 bytes, 5645 records. 0 BCF records skipped due to caller-specific exceptions
[1485789] [2023-09-08 12:31:10.919] [GLnexus] [info] Created sample set *@6
[1485789] [2023-09-08 12:31:10.919] [GLnexus] [info] Flushing database...
[1485789] [2023-09-08 12:31:11.545] [GLnexus] [info] Bulk load complete!
[1485789] [2023-09-08 12:31:11.558] [GLnexus] [info] found sample set *@6
[1485789] [2023-09-08 12:31:11.558] [GLnexus] [info] discovering alleles in 1 range(s) on 126 threads
[1485789] [2023-09-08 12:31:14.064] [GLnexus] [info] discovered 258742 alleles
[1485789] [2023-09-08 12:31:14.469] [GLnexus] [info] unified to 117841 sites cleanly with 122084 ALT alleles. 66 ALT alleles were additionally included in monoallelic sites and 8061 were filtered out on quality thresholds.
[1485789] [2023-09-08 12:31:14.469] [GLnexus] [info] Finishing database compaction...
[1485789] [2023-09-08 12:31:14.498] [GLnexus] [info] genotyping 117841 sites; sample set = *@6 mem_budget = 0 threads = 128
[1485789] [2023-09-08 12:31:20.343] [GLnexus] [info] genotyping complete!
[1485789] [2023-09-08 12:31:20.343] [GLnexus] [info] worker threads were cumulatively stalled for 456500ms
[1485789] [2023-09-08 12:31:20.343] [GLnexus] [info] Num BCF records read 4574092  query hits 727711

[user@cn3144 ~]$ **ls**
dv_platinum6_chr21.bcf	dv_platinum6_chr21_gvcf  GLnexus.DB  hg19_chr21.bed
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. glnexus.sh). For example:



```

#!/bin/bash
set -e
module load glnexus
cd /data/user
glnexus_cli --config DeepVariant --bed hg19_chr21.bed dv_platinum6_chr21_gvcf/*.gvcf.gz > dv_platinum6_chr21.bcf

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] glnexus.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. glnexus.swarm). For example:



```

glnexus_cli --config DeepVariant --bed genomic_range1.bed dv_platinum6_chr21_gvcf/*.gvcf.gz > output1.bcf
glnexus_cli --config DeepVariant --bed genomic_range2.bed dv_platinum6_chr21_gvcf/*.gvcf.gz > output2.bcf
glnexus_cli --config DeepVariant --bed genomic_range3.bed dv_platinum6_chr21_gvcf/*.gvcf.gz > output3.bcf
glnexus_cli --config DeepVariant --bed genomic_range4.bed dv_platinum6_chr21_gvcf/*.gvcf.gz > output4.bcf

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f glnexus.swarm [-g #] [-t #] --module glnexus
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module glnexus Loads the glnexus module for each subjob in the swarm 
 | |
 | |
 | |








