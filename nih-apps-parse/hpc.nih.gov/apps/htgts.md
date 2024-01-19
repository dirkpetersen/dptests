

document.querySelector('title').textContent = 'Htgts on HPC';
Htgts on HPC


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

  HTGTS represents High-Throughput Genome-Wide Translocation Sequencing pipeline 
 - provided by the [Alt 
 Lab](http://www.idi.harvard.edu/investigators_research/investigator/alt_lab/)

 ### 


Documentation * <http://robinmeyers.github.io/transloc_pipeline/thedocs.html>



Important Notes * Module Name: htgts (see [the modules 
 page](/apps/modules.html) for more information)
* Multithreaded





Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=20g --cpus-per-task=6**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load htgts**
[user@cn3144 ~]$ **GENOME\_DB=/fdb/htgts/genomes**
[user@cn3144 ~]$ **BOWTIE2\_INDEXES=/fdb/htgts/genomes/bowtie2\_indexes/hg19**
[user@cn3144 ~]$ **SPECIES=hg19**
[user@cn3144 ~]$ **cp -r /usr/local/apps/htgts/tutorial\_data/ /data/$USER/**
[user@cn3144 ~]$ **cd /data/$USER/tutorial\_data**
[user@cn3144 ~]$ **TranslocPreprocess.pl tutorial\_metadata.txt preprocess/ --threads $SLURM\_CPUS\_PER\_TASK --read1 pooled\_R1.fastq.gz --read2 pooled\_R2.fastq.gz**
[user@cn3144 ~]$ **TranslocWrapper.pl tutorial\_metadata.txt preprocess/ results/ --threads $SLURM\_CPUS\_PER\_TASK**
[user@cn3144 ~]$ **TranslocFilter.pl results/RAG1A\_IR\_SRep1/RAG1A\_IR\_SRep1.tlx results/RAG1A\_IR\_SRep1/RAG1A\_IR\_SRep1\_filtered.tlx --filters "f.unaligned=1"**
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
module load htgts
TranslocPreprocess.pl tutorial_metadata.txt preprocess/ --threads $SLURM_CPUS_PER_TASK --read1 poole\
d_R1.fastq.gz --read2 pooled_R2.fastq.gz                                                                                     TranslocWrapper.pl tutorial_metadata.txt preprocess/ results/ --threads $SLURM_CPUS_PER_TASK
TranslocFilter.pl results/RAG1A_IR_SRep1/RAG1A_IR_SRep1.tlx results/RAG1A_IR_SRep1/RAG1A_IR_SRep1_fi\
ltered.tlx --filters "f.unaligned=1"

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=6 --mem=20g batch.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. job.swarm). For example:



```

cd dir1; TranslocPreprocess.pl...;TranslocWrapper.pl...;TranslocFilter.pl...
cd dir2; TranslocPreprocess.pl...;TranslocWrapper.pl...;TranslocFilter.pl...
cd dir3; TranslocPreprocess.pl...;TranslocWrapper.pl...;TranslocFilter.pl...

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f job.swarm -g 20 -t 6 --module htgts
```

where
 

|  |  |
| --- | --- |
| -g *#*  | Number of Gigabytes of memory required for each process (1 line in the swarm command file)
  |
| -t *#* | Number of threads/CPUs required for each process (1 line in the swarm command file).
  |
| --module  | Loads the module for each subjob in the swarm  |




