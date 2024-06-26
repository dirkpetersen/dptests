

document.querySelector('title').textContent = 'pilon on Biowulf';
pilon on Biowulf


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


 Pilon is a application to improve draft assemblies and find variation among strains, including large event detection.



### References:


* Bruce J. Walker, Thomas Abeel, Terrance Shea, Margaret Priest, Amr Abouelliel, Sharadha Sakthikumar, Christina A. Cuomo, Qiandong Zeng, Jennifer Wortman, Sarah K. Young, Ashlee M. Earl.
*[Pilon: an integrated tool for comprehensive microbial variant detection and genome assembly improvement.](https://pubmed.ncbi.nlm.nih.gov/25409509/ )* . PLoS One. 2014 Nov 19;9(11).


Documentation
* pilon Main Site: [GitHub](https://github.com/broadinstitute/pilon)


Important Notes
* Module Name: pilon (see [the modules page](/apps/modules.html) for more information)
* Limit memory usage with the java heap size option -Xmx (e.g. -Xmx16g for 16 GB).
 * Environment variables set: *PILON\_HOME*, *PILON\_JAR** Note that to be able to use the environment variable *SLURM\_MEM\_PER\_NODE* below, you will need to explicitly allocate memory resources.



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --cpus-per-task=4 --mem=4g**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn4224 are ready for job

[user@cn4224 ~]$ **module load pilon**
[+] Loading pilon  1.23 

[user@cn4224 ~]$ **java -Xmx${SLURM\_MEM\_PER\_NODE}m -jar $PILON\_JAR --genome genome.fasta --frags frag.sorted.bam --ouput OUTPUT --vcf**
Pilon version 1.23 Mon Dec 21 10:00:01 2020 -0600
[...]

[user@cn4224 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226

[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. pilon.sh) similar to the following.



```

#! /bin/bash

set -e

module load pilon

java -Xmx${SLURM_MEM_PER_NODE}m -jar ${PILON_JAR} --genome genome.fasta --frags frag.sorted.bam --ouput OUTPUT --vcf

```

Submit these jobs using the Slurm [sbatch](/docs/userguide.html) command.


Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile for the first step of the pipeline (e.g. pilon.swarm). For example:



```

java -Xmx${SLURM_MEM_PER_NODE}m -jar ${PILON_JAR}  --genome genome.fasta --frags frag00.sorted.bam --ouput OUTPUT --vcf
java -Xmx${SLURM_MEM_PER_NODE}m -jar ${PILON_JAR}  --genome genome.fasta --frags frag02.sorted.bam --ouput OUTPUT --vcf
java -Xmx${SLURM_MEM_PER_NODE}m -jar ${PILON_JAR}  --genome genome.fasta --frags frag03.sorted.bam --ouput OUTPUT --vcf

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f pilon.swarm [-g #] --module pilon
```

where


|  |  |  |  |
| --- | --- | --- | --- |
| -g #  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | --module pilon  Loads the pilon module for each subjob in the swarm
 | |
 | |








