

document.querySelector('title').textContent = 'EukDetect on Biowulf';
EukDetect on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
 |



EukDetect is a python-based workflow application that detects eukaryotes in shotgun metagenomic data using eukaryotic gene markers. According to the authors: 




> 
>  EukDetect is accurate, sensitive, has a broad taxonomic coverage of microbial eukaryotes, and is resilient against bacterial contamination in eukaryotic genomes.
> 



The application can be run using the eukdetect executable or with snakemake that comes bundled with the application.
The gene markers are available in /fdb/eukdetect/eukdb. User must create a configuration YAML for the pipeline.
A default YAML is available to Biowulf users with certain fields edited based on the installation.



### References:


* Lind, A. L. and K.S. Pollard.
 [**Accurate and sensitive detection of microbial eukaryotes from metagenomic shotgun sequencing data.**](https://doi.org/10.1101/2020.07.22.216580)
*bioRxiv. 2020.07.22.216580. (2020).*


Documentation
* [EukDetect GitHub](https://github.com/allind/EukDetect)


Important Notes
* Module Name: eukdetect (see [the modules page](/apps/modules.html) for more information)
* Multithreaded
* Environment variables set 
	+ EUKDETECT\_DB
	+ EUKDETECT\_EXAMPLES
	+ EUKDETECT\_SHARE* Example fastq in /usr/local/apps/eukdetect/examples
* Reference data in /fdb/eukdetect/eukdb
* Default configuration YAML in /usr/local/apps/eukdetect/1.2/share



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive --gres=lscratch:10 --mem=8G -c4**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load eukdetect**

[user@cn3144 ~]$ **cd /lscratch/$SLURM\_JOB\_ID**

[user@cn3144 ~]$ **mkdir fastq**

[user@cn3144 ~]$ **cp $EUKDETECT\_EXAMPLES/\*fastq.gz fastq/**

[user@cn3144 ~]$ **cp $EUKDETECT\_SHARE/default\_configfile.yml .**

[user@cn3144 ~]$ **nano default\_configfile.yml**
# Edit YAML and set output_dir, fq_dir, and samples fields correctly

[user@cn3144 ~]$ **eukdetect --mode runall --configfile default\_configfile.yml --cores 4**
  01/20/2022 16:43:16:  Parsing config file ...
  01/20/2022 16:43:16:  Running: snakemake --snakefile /opt/EukDetect/rules/eukdetect.rules --configfile default_configfile.yml --cores 4
  01/20/2022 16:43:16:  Redirecting snakemake output to snakemake_1642714996.749603.log
  01/20/2022 16:43:46:  Snakemake finished
  01/20/2022 16:43:46:  Snakemake pipeline created all files.

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. eukdetect.sh). For example:



```

#!/bin/bash
set -e
module load eukdetect
eukdetect --mode runall --configfile default_configfile.yml --cores $SLURM_CPUS_PER_TASK

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] eukdetect.sh
```







