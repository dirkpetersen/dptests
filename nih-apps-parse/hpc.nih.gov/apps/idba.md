

document.querySelector('title').textContent = 'IDBA on Biowulf';
IDBA on Biowulf


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



IDBA is the basic iterative de Bruijn graph assembler for second-generation sequencing reads. IDBA-UD, an extension of IDBA, is designed to utilize paired-end reads to assemble low-depth regions and use progressive depth on contigs to reduce errors in high-depth regions. It is a generic purpose assembler and especially good for single-cell and metagenomic sequencing data. IDBA-Hybrid is another update version of IDBA-UD, which can make use of a similar reference genome to improve assembly result. IDBA-Tran is an iterative de Bruijn graph assembler for RNA-Seq data.


The basic IDBA is included only for comparison.
* If you are assembling genomic data without reference, please use IDBA-UD.
* If you are assembling genomic data with a similar reference genome, please use IDBA-Hybrid.
* If you are assembling transcriptome data, please use IDBA-Tran.





### References:


* Peng, Y., et al. (2010) IDBA- A Practical Iterative de Bruijn Graph De Novo Assembler. RECOMB. Lisbon. doi: [10.1007/978-3-642-12683-3\_28](https://doi.org/10.1007/978-3-642-12683-3_28)* Peng, Y., et al. (2012) IDBA-UD: a de novo assembler for single-cell and metagenomic sequencing data with highly uneven depth, Bioinformatics, 28, 1420-1428.Paper, doi: [10.1093/bioinformatics/bts174](https://doi.org/10.1093/bioinformatics/bts174)


Documentation
* [IDBA suite on GitHub](https://github.com/loneknightpy/idba)* [IDBA](http://i.cs.hku.hk/~alse/hkubrg/projects/idba/index.html)* [IDBA-UD](http://i.cs.hku.hk/~alse/hkubrg/projects/idba_ud/index.html)* [IDBA-Hybrid](http://i.cs.hku.hk/~alse/hkubrg/projects/idba_hybrid/index.html)* [IDBA-Tran](http://i.cs.hku.hk/~alse/hkubrg/projects/idba_tran/index.html)


Important Notes
* Module Name: idba (see [the modules page](/apps/modules.html) for more information)

Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job
[user@cn3144 ~]$ **module load idba**
[user@cn3144 ~]$ **idba\_ud -r reads.fa -o output**
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. TEMPLATE.sh). For example:



```

#!/bin/bash
set -e
module load idba
idba_ud  -r read.fa -o output

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] TEMPLATE.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. idba.swarm). For example:



```

idba_ud  -r sample1-reads.fa -o sample1
idba_ud  -r sample2-reads.fa -o sample2
idba_ud  -r sample3-reads.fa -o sample3
idba_ud  -r sample4-reads.fa -o sample4

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f idba.swarm [-g #] [-t #] --module idba
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module idba Loads the IDBA module for each subjob in the swarm 
 | |
 | |
 | |


