

document.querySelector('title').textContent = 'Hisat on Biowulf';
Hisat on Biowulf


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



HISAT is a fast and sensitive spliced alignment program. As part of HISAT, it includes a new indexing scheme based on the Burrows-Wheeler transform (BWT) and the FM index, called hierarchical indexing, that employs two types of indexes: (1) one global FM index representing the whole genome, and (2) many separate local FM indexes for small regions collectively covering the genome. ned for mapping different types of RNA-seq reads. All these together, HISAT enables extremely fast and sensitive alignment of reads, in particular those spanning two exons or more. HISAT is distributed under theGPLv3 license.



Beside the index files under **/fdb/hisat**, HISAT on Biowulf has been built with sra support. A local copy of the SRA Refseq reference data is maintained in **/fdb/sra/refseq**, and 
the NCBI SRA toolkit has been configured on Biowulf to :
1. first check if the required reference data is available in the local copy
- if not, download the required data using the Biowulf web proxy


Thus, running hisat on the Biowulf compute nodes with SRA reference data should be seamless. 

### References:


* Kim et al. HISAT: a fast spliced aligner with low memory requirements. Nature Methods 12, 357-360 (2015)


Documentation
* [Hisat Main Site](http://www.ccb.jhu.edu/software/hisat/manual.shtml)



Important Notes
* Module Names: The module naming system for hisat can be confusing, so the table below is provided for clarity:


|  |  |  |  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| **Version  **Module Name  **Alternate Module Name 
| Hisat2 v 2.2.1 built for NCBI NGS v3.0.1 hisat/2.2.2.1-ngs3.0.1 hisat2/2.2.1 
| Hisat2 v 2.1.0 built for NCBI NGS v2.10.9 hisat/2.2.1.0-ngs2.10.9  hisat2/2.1.0 

 | | |
 | | |** |** |** |


* Multithreaded app
* environment variables set 
	+ $HISAT\_INDEXES
	+ $HISAT\_HOME* Example files in $HISAT\_HOME/example* Reference data in /fdb/hisat/



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --cpus-per-task=4 --mem=10g**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load hisat**

[user@cn3144 ~]$ **cd /data/$USER/hisat**

[user@cn3144 ~]$ **hisat2 -p 4 -x $HISAT\_INDEXES/grch38/genome -f $HISAT\_HOME/example/reads/reads\_1.fa \**
**-S /data/$USER/output.sam** 

[user@cn3144 ~]$ **hisat2 -p 4 -x $HISAT\_INDEXES/grch38/genome -f $HISAT\_HOME/example/reads/reads\_1.fa \**
**$HISAT\_HOME/example/reads/reads\_2.fa -S /data/$USER/output.sam**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. hisat.sh). For example:



```

#!/bin/bash
cd /data/$USER/mydir
module load hisat

hisat2 -p $SLURM_CPUS_PER_TASK -x $HISAT_INDEXES/grch38/genome -f $HISAT_HOME/example/reads/reads_1.fa \
-S /data/$USER/output1.sam
hisat2 -p $SLURM_CPUS_PER_TASK -x $HISAT_INDEXES/grch38/genome -f $HISAT_HOME/example/reads/reads_1.fa \
 $HISAT_HOME/example/reads/reads_2.fa -S /data/$USER/output2.sam
......
......

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=4 --mem=10g hisat.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. hisat.swarm). For example:



```

cd dir1; hisat2 -p $SLURM_CPUS_PER_TASK -x $HISAT_INDEXES/grch38/genome \
-f $HISAT_HOME/example/reads/reads_1.fa -S /data/$USER/output1.sam
cd dir2; hisat2 -p $SLURM_CPUS_PER_TASK -x $HISAT_INDEXES/grch38/genome \
-f $HISAT_HOME/example/reads/reads_1.fa -S /data/$USER/output1.sam
cd dir3; hisat2 -p $SLURM_CPUS_PER_TASK -x $HISAT_INDEXES/grch38/genome \
-f $HISAT_HOME/example/reads/reads_1.fa -S /data/$USER/output1.sam
[...]   

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f hisat.swarm [-g 10] [-t 4] --module hisat
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module hisat Loads the hisat module for each subjob in the swarm 
 | |
 | |
 | |


Benchmarks
This benchmarks used 3gb paired sequences. 16 cores is most efficient.


 




| Core | Time(hr:min) |
| --- | --- |
| 1 | 1:45:00 |
| 2 | 1:05:00 |
| 4 | 39:50 |
| 8 | 21:51 |
| 16 | 12:46 |
| 32 | 11:39 |















