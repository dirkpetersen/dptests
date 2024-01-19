

document.querySelector('title').textContent = 'Bamliquidator on Biowulf';
Bamliquidator on Biowulf


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



Bamliquidator is a set of tools for analyzing the density of short DNA sequence read alignments in the BAM file format

* the read counts across multiple genomes are grouped, normalized, summarized, and graphed in interactive html files
* for an interactive graph example, see this summary and this breakdown for a single chromosome
* a BAM file is a binary sequence alignment map -- see SAMtools for more info
* the read counts and summaries are stored in HDF5 format where they can be efficiently read via Python PyTables or the HDF5 C apis
+ the HDF5 files can be viewed directly with the cross platform tool HDFView
+ there is an option to output the data in tab delimited text files as well

* there is also a command line utility for counting the number of reads in specified portion of a chromosome, and the count is output to the console





Documentation
* [bamliquidator Main Site](https://github.com/BradnerLab/pipeline/wiki/bamliquidator)


Important Notes
* Module Name: bamliquidator (see [the modules page](/apps/modules.html) for more information)
* Multithreaded app



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

[user@cn3144 ~]$ **cd /lscratch/$SLURM\_JOB\_ID**

[user@cn3144 ~]$ **cp -r /usr/local/apps/bamliquidator/test .**

[user@cn3144 ~]$ **module load bamliqudiator**

[user@cn3144 ~]$ **bamliquidator\_batch test/wgEncodeUwRepliSeqBg02esG1bAlnRep1.bam**
Liquidating test/wgEncodeUwRepliSeqBg02esG1bAlnRep1.bam (file 1 of 1)
Liquidation completed: 4.188048 seconds, 1784867 reads, 0.238775 millions of reads per second
Cell Types: test
Normalizing and calculating percentiles for cell type test
Indexing normalized counts
Plotting
-- skipping plotting chrM because not enough bins (only 1)
-- skipping plotting chrM because not enough bins (only 1)
Summarizing
Post liquidation processing took 2.364760 seconds

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. bamliquidator.sh). For example:



```

#!/bin/bash
module load bamliquidator
bamliquidator_batch.py -n 16 -o test.output directory_of_bam_files/

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=16 --mem=10g bamliquidator.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. bamliquidator.swarm). For example:



```

bamliquidator sorted_and_indexed.bam chr1 0 1000000 + 2 200
bamliquidator sorted_and_indexed.bam chr1 1000000 2000000 + 2 200
bamliquidator sorted_and_indexed.bam chr1 2000000 3000000 + 2 200
bamliquidator sorted_and_indexed.bam chr1 3000000 4000000 + 2 200

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f bamliquidator.swarm [-g 4] [-t 2] --module bamliquidator
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module bamliquidator Loads the bamliquidator module for each subjob in the swarm 
 | |
 | |
 | |








