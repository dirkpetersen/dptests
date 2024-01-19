

document.querySelector('title').textContent = 'Pyclone on HPC';
Pyclone on HPC


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

  PyClone is statistical model and software tool designed to infer the prevalence 
 of point mutations in heterogeneous cancer samples. The input data for PyClone 
 consists of a set read counts from a deep sequencing experiment, the copy 
 number of the genomic region containing the mutation and an estimate of 
 tumour content. 


### References:


* Roth *et al.* PyClone: statistical inference of clonal population 
 structure in cancer [PMID: 
 24633410](http://www.nature.com/nmeth/journal/v11/n4/full/nmeth.2883.html)


Documentation * <https://github.com/aroth85/pyclone>



Important Notes * Module Name: pyclone (see [the 
 modules page](/apps/modules.html) for more information)
* Example files in /usr/local/apps/pyclone/examples





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

[user@cn3144 ~]$ **module load pyclone**
[user@cn3144 ~]$ **cp -r /usr/local/apps/pyclone/examples /data/$USER**
[user@cn3144 ~]$ **cd /data/$USER/examples/mixing/tsv**
[user@cn3144 ~]$ **PyClone run\_analysis\_pipeline --in\_files SRR385939.tsv SRR385940.tsv SRR385941.tsv --working\_dir pyclone\_analysis**
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```

This will create a directory pyclone\_analysis. After the command completes the directory will contain several folders and the file config.yaml:
* config.yaml
* plots/
* tables/
* trace/
* yaml/


The contents of these folders are as follows:
* config.yaml - This file specifies the configuration used for the PyClone analysis.
* plots - Contains all plots from the analysis. There will be two sub-folders clusters/ and loci/ for cluster and locus specific plots respectively.
* tables - This contains the output tables with summarized results for the analysis. There will be two tables clusters.tsv and loci.tsv, for cluster and locus specific information.
* trace - This the raw trace from the MCMC sampling algorithm. Advanced users may wish to work with these files directly for generating plots and summary statistics.




Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. batch.sh). For example:



```

#!/bin/bash
set -e
module load pyclone
PyClone run_analysis_pipeline --in_files SRR385939.tsv SRR385940.tsv SRR385941.tsv --working_dir pyclone_analysis
```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--mem=#] batch.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. job.swarm). For example:



```

cd dir1;PyClone ...
cd dir2;PyClone ...
cd dir3;PyClone ...
cd dir4;PyClone ...

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f job.swarm [-g #] --module pyclone
```

where
 

|  |  |
| --- | --- |
| -g *#*  | Number of Gigabytes of memory required for each process 
 (1 line in the swarm command file)  |
| --module  | Loads the module for each subjob in the swarm  |








