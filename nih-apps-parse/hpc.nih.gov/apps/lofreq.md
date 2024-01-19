

document.querySelector('title').textContent = 'Lofreq on HPC';
Lofreq on HPC


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

  LoFreq (i.e. LoFreq version 2) is a fast and sensitive variant-caller for 
 inferring SNVs and indels from next-generation sequencing data. It makes 
 full use of base-call qualities and other sources of errors inherent in 
 sequencing (e.g. mapping or base/indel alignment uncertainty), which are 
 usually ignored by other methods or only used for filtering. 


LoFreq can run on almost any type of aligned sequencing data (e.g. Illumina, 
 IonTorrent or Pacbio) since no machine- or sequencing-technology dependent 
 thresholds are used. It automatically adapts to changes in coverage and 
 sequencing quality and can therefore be applied to a variety of data-sets 
 e.g. viral/quasispecies, bacterial, metagenomics or somatic data.


LoFreq is very sensitive; most notably, it is able to predict variants 
 below the average base-call quality (i.e. sequencing error rate). Each variant 
 call is assigned a p-value which allows for rigorous false positive control. 
 Even though it uses no approximations or heuristics, it is very efficient 
 due to several runtime optimizations and also provides a (pseudo-)parallel 
 implementation. LoFreq\* is generic and fast enough to be applied to high-coverage 
 data and large genomes. On a single processor it takes a minute to analyze 
 Dengue genome sequencing data with nearly 4000X coverage, roughly one hour 
 to call SNVs on a 600X coverage E.coli genome and also roughly an hour to 
 run on a 100X coverage human exome dataset.


### References:

 * <http://www.ncbi.nlm.nih.gov/pubmed/23066108>


Documentation * <http://csb5.github.io/lofreq/commands/>



Important Notes * Module Name: lofreq (see [the modules 
 page](/apps/modules.html) for more information)





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

[user@cn3144 ~]$ **module load lofreq**
[user@cn3144 ~]$ **lofreq call -f ref.fa -o vars.vcf aln.bam**

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
module load lofreq
lofreq call -f ref.fa -o vars.vcf aln.bam
```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--mem=#] batch.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. job.swarm). For example:



```
cd dir1; lofreq call -f ref.fa -o vars.vcf aln.bam
cd dir2; lofreq call -f ref.fa -o vars.vcf aln.bam
cd dir3; lofreq call -f ref.fa -o vars.vcf aln.bam
```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f job.swarm [-g #] --module lofreq
```

where
 

|  |  |
| --- | --- |
| -g *#*  | Number of Gigabytes of memory required for each process 
 (1 line in the swarm command file)  |
| --module  | Loads the module for each subjob in the swarm  |




