

document.querySelector('title').textContent = 'pbipa on Biowulf';
pbipa on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
 |



pbipa is a genome assembler for PacBio HiFi reads.

> 
>  Improved Phased Assembler (IPA) is the official PacBio software for
>  HiFi genome assembly. IPA was designed to utilize the accuracy of 
>  PacBio HiFi reads to produce high-quality phased genome assemblies.
>  IPA is an end-to-end solution, starting with input reads and resulting
>  in a polished assembly. IPA is fast, providing an easy to use local
>  run mode or a distributed pipeline for a cluster.
> 


Under the hood, pbipa uses snakemake to accomplish the various assembly
steps from phasing and filtering to contig construction and read polishing.

The installation on Biowulf also provides the following applications:
* [falconc](https://bioconda.github.io/recipes/pb-falconc/README.html)
* [racon](https://github.com/lbcb-sci/racon)
* [nighthawk](https://www.pacb.com/blog/direct-phased-genome-assembly-using-nighthawk-on-hifi-reads/)


Documentation
* [GitHub Wiki](https://github.com/PacificBiosciences/pbbioconda/wiki/Improved-Phased-Assembler)


Important Notes
* Module Name: pbipa (see [the modules page](/apps/modules.html) for more information)
* Please refrain from loading snakemake and other python modules when running pbipa
* Multithreaded


We do not currently support pbipa to be run in distributed (dist) mode. Please only run in local mode.



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run ipa. Here we show how to run in IPA's local mode.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive --cpus-per-task=8 --mem=12G --gres=lscratch:20**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load pbipa**

[user@cn3144 ~]$ **ipa local --nthreads 8 --njobs 1 -i hifi.fasta.gz**

INFO: /opt/conda/envs/ipa/bin/ipa local --nthreads 8 --njobs 1 -i hifi.fasta.gz
INFO: ipa.py ipa (wrapper) version=1.5.0 ... Checking dependencies ...
INFO: Dependencies
/opt/conda/envs/ipa/bin/python3
/opt/conda/envs/ipa/bin/ipa2-task
/opt/conda/envs/ipa/bin/falconc
/opt/conda/envs/ipa/bin/minimap2
/opt/conda/envs/ipa/bin/nighthawk
/opt/conda/envs/ipa/bin/pancake
/opt/conda/envs/ipa/bin/pblayout
/opt/conda/envs/ipa/bin/racon
/opt/conda/envs/ipa/bin/samtools
/opt/conda/envs/ipa/bin/ipa_purge_dups
/opt/conda/envs/ipa/bin/ipa_purge_dups_split_fa
snakemake version=6.8.1
ipa2-task 1.5.0 (commit c875fce13bdacbafc2f4f750c6438f4453e1354d)
 Machine name: 'Linux'
Copyright (C) 2004-2021     Pacific Biosciences of California, Inc.
This program comes with ABSOLUTELY NO WARRANTY; it is intended for
Research Use Only and not for use in diagnostic procedures.

falconc version=1.13.1+git.f9d1b5651e891efe379bd9727a0fa0931b875d7b, nim-version=1.5.1
minimap2 version=2.22-r1101
Nighthawk 0.1.0 (commit SL-release-10.1.0-7-gbe5dfb1*)
pancake 1.3.0 (commit SEQII-release-10.1.0-432-gf2693fd*)
pblayout 1.0.0 (commit SL-release-10.1.0-152-g66936d1*)
racon version=v1.4.20
samtools 1.12
Using htslib 1.12
ipa_purge_dups Version: 1.2.5
...

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. ipa.sh). Again, this shows how to setup a job in 'local' mode. For example:



```

#!/bin/bash
set -e
module load pbipa
ipa local --nthreads $SLURM_CPUS_PER_TASK --njobs 1 -i input.fasta

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] ipa.sh
```













