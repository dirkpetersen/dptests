

document.querySelector('title').textContent = 'Magic-BLAST on Biowulf';
Magic-BLAST on Biowulf


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



Magic-BLAST is a tool for mapping large next-generation RNA or DNA sequencing runs against a whole genome or transcriptome. Each alignment optimizes a composite score, taking into account simultaneously the two reads of a pair, and in case of RNA-seq, locating the candidate introns and adding up the score of all exons. This is very different from other versions of BLAST, where each exon is scored as a separate hit and read-pairing is ignored.




Magic-BLAST incorporates within the NCBI BLAST code framework ideas developed in the NCBI Magic pipeline, in particular [hit extensions by local walk and jump](http://www.ncbi.nlm.nih.gov/pubmed/26109056), and recursive clipping of mismatches near the edges of the reads, which avoids accumulating artefactual mismatches near splice sites and is needed to distinguish short indels from substitutions near the edges.



### References:


* Boratyn GM, Thierry-Mieg J, Thierry-Mieg D, Busby B, Madden TL. (2019) **Magic-BLAST, an accurate RNA-seq aligner for long and short reads.** BMC Bioinformatics 20: 405. [doi:10.1186/s12859-019-2996-x](https://doi.org/10.1186/s12859-019-2996-x)


Documentation
* [Magic-BLAST Main Site](https://ncbi.github.io/magicblast)


Important Notes
* Module Name: magicblast (see [the modules page](/apps/modules.html) for more information)
 * Multithreaded using the -num\_threads argument to magicblast.



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.
  
The following sample session (user input in **bold**) is based on cookbooks from [the official documentation](https://ncbi.github.io/magicblast/):



```

[user@biowulf]$ **sinteractive --cpus-per-task 8**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load magicblast**

[user@cn3144 ~]$ **makeblastdb -in my\_reference.fa -out my\_reference -parse\_seqids -dbtype nucl** 
[user@cn3144 ~]$ **magicblast -num\_threads $SLURM\_CPUS\_PER\_TASK -query reads.fa -db my\_reference**
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. magicblast.sh). For example:



```

#!/bin/bash
set -e
module load magicblast
test -n $SLURM_CPUS_PER_TASK || SLURM_CPUS_PER_TASK=2
magicblast -num_threads $SLURM_CPUS_PER_TASK -query reads.fa -db my_reference

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] magicblast.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. magicblast.swarm). For example:



```

cd sample1 && magicblast -num_threads $SLURM_CPUS_PER_TASK -query reads.fa -db /data/$USER/references/my_reference
cd sample2 && magicblast -num_threads $SLURM_CPUS_PER_TASK -query reads.fa -db /data/$USER/references/my_reference
cd sample3 && magicblast -num_threads $SLURM_CPUS_PER_TASK -query reads.fa -db /data/$USER/references/my_reference
cd sample4 && magicblast -num_threads $SLURM_CPUS_PER_TASK -query reads.fa -db /data/$USER/references/my_reference

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f magicblast.swarm [-g #] -t # --module magicblast
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module magicblast Loads the Magic-BLAST module for each subjob in the swarm 
 | |
 | |
 | |








