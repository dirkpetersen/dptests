

document.querySelector('title').textContent = 'Somaticsniper on HPC';
Somaticsniper on HPC


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

  The purpose of this program is to identify single nucleotide positions that 
 are different between tumor and normal (or, in theory, any two bam files). 
 It takes a tumor bam and a normal bam and compares the two to determine 
 the differences. It outputs a file in a format very similar to Samtools 
 consensus format. It uses the genotype likelihood model of MAQ (as implemented 
 in Samtools) and then calculates the probability that the tumor and normal 
 genotypes are different. This probability is reported as a somatic score. 
 The somatic score is the Phred-scaled probability (between 0 to 255) that 
 the Tumor and Normal genotypes are not different where 0 means there is 
 no probability that the genotypes are different and 255 means there is a 
 probability of 1 – 10(255/-10) that the genotypes are different between 
 tumor and normal. This is consistent with how the SAM format reports such 
 probabilities. It is currently available as source code via github or as 
 a Debian APT package. 


This tool is developed by [David 
 E. Larson etc](http://gmt.genome.wustl.edu/somatic-sniper/current/).


### References:

 * <https://academic.oup.com/bioinformatics/article/28/3/311/188933>


Documentation * <http://gmt.genome.wustl.edu/packages/somatic-sniper/>



Important Notes * Module Name: somaticsniper (see [the 
 modules page](/apps/modules.html) for more information)





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

[user@cn3144 ~]$ **module load somaticsniper**
[user@cn3144 ~]$ **bam-somaticsniper -f ref.fasta tumor.bam normal.bam Outfile**

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
module load somaticsniper
bam-somaticsniper -f ref.fasta tumor.bam normal.bam Outfile
```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--mem=#] batch.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. job.swarm). For example:



```

cd dir1; bam-somaticsniper -f ref.fasta tumor.bam normal.bam Outfile
cd dir2; bam-somaticsniper -f ref.fasta tumor.bam normal.bam Outfile
cd dir3; bam-somaticsniper -f ref.fasta tumor.bam normal.bam Outfile
cd dir4; bam-somaticsniper -f ref.fasta tumor.bam normal.bam Outfile


```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f job.swarm [-g #] --module somaticsniper
```

where
 

|  |  |
| --- | --- |
| -g *#*  | Number of Gigabytes of memory required for each process 
 (1 line in the swarm command file)  |
| --module  | Loads the module for each subjob in the swarm  |




