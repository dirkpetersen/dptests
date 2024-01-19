

document.querySelector('title').textContent = 'Umitools on HPC';
Umitools on HPC


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

  Umi-tools are tools for dealing with Unique Molecular Identifiers (UMIs)/Random 
 Molecular Tags (RMTs) and single cell RNA-Seq cell barcodes. Currently there 
 are 6 commands.


The extract and whitelist commands are used to prepare a fastq containg 
 UMIs +/- cell barcodes for alignment.


whitelist:  

 Builds a whitelist of the 'real' cell barcodes  

 This is useful for droplet-based single cell RNA-Seq where the identity 
 of the true cell barcodes is unknown. Whitelist can then be used to filter 
 with extract (see below)


extract:  

 Flexible removal of UMI sequences from fastq reads.  

 UMIs are removed and appended to the read name. Any other barcode, for example 
 a library barcode, is left on the read. Can also filter reads by quality 
 or against a whitelist (see above)  

 The remaining commands, group, dedup and count/count\_tab, are used to identify 
 PCR duplicates using the UMIs and perform different levels of analysis depending 
 on the needs of the user. A number of different UMI deduplication schemes 
 are enabled - The recommended method is directional.


group:  

 Groups PCR duplicates using the same methods available through `dedup`.  

 This is useful when you want to manually interrogate the PCR duplicates


dedup:  

 Groups PCR duplicates and deduplicates reads to yield one read per group  

 Use this when you want to remove the PCR duplicates prior to any downstream 
 analysis


count:  

 Groups and deduplicates PCR duplicates and counts the unique molecules per 
 gene  

 Use this when you want to obtain a matrix with unique molecules per gene. 
 Can also perform per-cell counting for scRNA-Seq.


count\_tab:  

 As per count except input is a flatfile


### 


Documentation * <https://github.com/CGATOxford/UMI-tools>



Important Notes
* Module Name: umitools (see [the 
 modules page](/apps/modules.html) for more information)
* Example files in /usr/local/apps/umitools/example.fastq.gz





Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load umitools**
[user@cn3144 ~]$ **cp /usr/local/apps/umitools/example.fastq.gz .** 
[user@cn3144 ~]$ **umi\_tools extract --stdin=example.fastq.gz --bc-pattern=NNNNNNNNN --log=processed.log --stdout processed.fastq.gz**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```




Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. batch.sh). For example:



```

#!/bin/bash
set -e
module load umitools
umi_tools extract --stdin=example.fastq.gz --bc-pattern=NNNNNNNNN --log=processed.log --stdout processed.fastq.gz
```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--mem=#] batch.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. job.swarm). For example:



```

cd dir1; umi_tools extract --stdin=example.fastq.gz --bc-pattern=NNNNNNNNN --log=processed.log --stdout processed.fastq.gz
cd dir2; umi_tools extract --stdin=example.fastq.gz --bc-pattern=NNNNNNNNN --log=processed.log --stdout processed.fastq.gz
cd dir3; umi_tools extract --stdin=example.fastq.gz --bc-pattern=NNNNNNNNN --log=processed.log --stdout processed.fastq.gz

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f job.swarm [-g #] --module umitools
```

where
 

|  |  |
| --- | --- |
| -g *#*  | Number of Gigabytes of memory required for each process 
 (1 line in the swarm command file)  |
| --module  | Loads the module for each subjob in the swarm  |




