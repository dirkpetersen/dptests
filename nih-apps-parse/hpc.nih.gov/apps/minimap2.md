

document.querySelector('title').textContent = 'minimap2 on Biowulf';
minimap2 on Biowulf


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



Minimap2 is a fast sequence mapping and alignment program that can find overlaps between long noisy reads, or map long reads or their assemblies to a reference genome optionally with detailed alignment (i.e. CIGAR). At present, it works efficiently with query sequences from a few kilobases to ~100 megabases in length at an error rate ~15%. Minimap2 outputs in the PAF or the SAM format. On limited test data sets, minimap2 is over 20 times faster than most other long-read aligners. It will replace BWA-MEM for long reads and contig alignment.



### References:


* Heng Li; Minimap2: pairwise alignment for nucleotide sequences, Bioinformatics, , bty191, doi:[10.1093/bioinformatics/bty191](https://doi.org/10.1093/bioinformatics/bty191)


Documentation
* [minimap2 main site](https://github.com/lh3/minimap2)
* minimap2(1)


Important Notes
* Module Name: minimap2 (see [the modules page](/apps/modules.html) for more information)
* environment variables set
	+ MINIMAP2\_HOME* Example files in $MINIMAP2\_HOME/test and /fdb/minimap2



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
[user@cn3144 ~]$ **module load minimap2**
[+] Loading minimap2, version 2.2...
[user@cn3144 ~]$ **man minimap2** # read the documentation
[user@cn3144 ~]$ **minimap2 -ax map10k $MINIMAP2\_HOME/test/MT-human.fa $MINIMAP2\_HOME/test/MT-orang.fa > test.sam**
[M::mm_idx_gen::0.005*1.16] collected minimizers
[M::mm_idx_gen::0.007*1.24] sorted minimizers
[M::main::0.007*1.24] loaded/built the index for 1 target sequence(s)
[M::mm_mapopt_update::0.008*1.19] mid_occ = 3
[M::mm_idx_stat] kmer size: 19; skip: 10; is_HPC: 1; #seq: 1
[M::mm_idx_stat::0.008*1.27] distinct minimizers: 2127 (99.95% are singletons); average occurrences: 1.000; average spacing: 7.786
[M::worker_pipeline::0.035*1.03] mapped 1 sequences
[M::main] Version: 2.2-r409
[M::main] CMD: minimap2 -ax map10k /usr/local/apps/minimap2/2.2/test/MT-human.fa /usr/local/apps/minimap2/2.2/test/MT-orang.fa
[M::main] Real time: 0.036 sec; CPU: 0.037 sec
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. minimap2.sh). For example:



```

#!/bin/bash
set -e
module load minimap2
minimap2 -ax map10k $MINIMAP2_HOME/test/MT-human.fa $MINIMAP2_HOME/test/MT-orang.fa > test.sam

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] minimap2.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. minimap2.swarm). For example:



```

minimap2 -ax map10k reference.fa reads1.fa > out1.sam
minimap2 -ax map10k reference.fa reads2.fa > out2.sam
minimap2 -ax map10k reference.fa reads3.fa > out3.sam
minimap2 -ax map10k reference.fa reads4.fa > out4.sam

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f minimap2.swarm [-g #] [-t #] --module minimap2
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g #  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t #  Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module minimap2  Loads the minimap2 module for each subjob in the swarm
 | |
 | |
 | |








