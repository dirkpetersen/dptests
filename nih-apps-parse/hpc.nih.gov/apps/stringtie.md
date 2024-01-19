

document.querySelector('title').textContent = 'Stringtie on Biowulf';
Stringtie on Biowulf


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


**StringTie** is a fast and highly efficient assembler of RNA-Seq alignments into potential transcripts. It uses a novel network flow algorithm as well as an optional *de novo* assembly step to assemble and quantitate full-length transcripts representing multiple splice variants for each gene locus. Its input can include not only the alignments of raw reads used by other transcript assemblers, but also alignments longer sequences that have been assembled from those reads.To identify differentially expressed genes between experiments, StringTie's output can be processed either by the [Cuffdiff](http://cole-trapnell-lab.github.io/cufflinks/cuffdiff/index.html) or [Ballgown](https://github.com/alyssafrazee/ballgown) programs.


StringTie is free, open source software released under an [Artistic Licen](http://opensource.org/licenses/artistic-license-2.0)


### References:

 * Pertea M, Pertea GM, Antonescu CM, Chang TC, Mendell JT & Salzberg 
 SL. [StringTie enables 
 improved reconstruction of a transcriptome from RNA-seq reads](https://www.nature.com/articles/nbt.3122) Nature 
 Biotechnology 2015, doi:10.1038/nbt.3122
* Pertea M, Kim D, Pertea GM, Leek JT, Salzberg SL [Transcript-level 
 expression analysis of RNA-seq experiments with HISAT, StringTie and Ballgown](https://www.nature.com/articles/nprot.2016.095), 
 Nature Protocols 11, 1650-1667 (2016), doi:10.1038/nprot.2016.095


Documentation * <https://ccb.jhu.edu/software/stringtie/index.shtml>



Important Notes * Module Name: stringtie (see [the 
 modules page](/apps/modules.html) for more information)
* Multithreaded





Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --cpus-per-task=4**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load stringtie**
[user@cn3144 ~]$ **stringtie <aligned\_reads.bam> -p $SLURM\_CPUS\_PER\_TASK [other options]**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```




Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit). Create a batch input file (e.g. stringtie.sh). For example:



```

#!/bin/bash
set -e
module load stringtie
stringtie <aligned_reads.bam> -p $SLURM_CPUS_PER_TASK [other options]
```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=4 --mem=10g stringtie.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. stringtie.swarm). For example:



```

cd dir1; stringtie aligned_reads.bam -p $SLURM_CPUS_PER_TASK [other options]
cd dir2; stringtie aligned_reads.bam -p $SLURM_CPUS_PER_TASK [other options]
cd dir3; stringtie aligned_reads.bam -p $SLURM_CPUS_PER_TASK [other options]

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f stringtie.swarm g 10 -t 4 --module stringtie
```

where
 

|  |  |
| --- | --- |
| -g *#*  | Number of Gigabytes of memory required for each process (1 line in the swarm command file)
  |
| -t *#* | Number of threads/CPUs required for each process (1 line in the swarm command file).
  |
| --module  | Loads the module for each subjob in the swarm  |




