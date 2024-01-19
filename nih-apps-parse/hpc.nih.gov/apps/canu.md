

document.querySelector('title').textContent = 'canu on Biowulf';
canu on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
[Batch job using Grid options](#grid) 
 |



Canu is a fork of the Celera Assembler designed for high-noise single-molecule sequencing (such as the PacBio RSII or Oxford Nanopore MinION). Canu will correct the reads, then trim suspicious regions (such as remaining SMRTbell adapter), then assemble the corrected and cleaned reads into unitigs. 



### References:


* Canu was developed by Adam Phillippy, Sergey Koren, Brian Walenz. 
* ([Canu website](http://canu.readthedocs.org/en/stable/))


Documentation
* [canu website](http://canu.readthedocs.org/en/stable)


Important Notes
* Module Name: canu (see [the modules page](/apps/modules.html) for more information)
* Multithreaded
* Test data for Canu is available in  
/usr/local/apps/canu/p6.25x.fastq (223 MB)



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive --cpus-per-task=8 --mem=20g**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load canu**
[+] Loading gnuplot 5.2.2  ...
[+] Loading canu  1.7

[user@cn3144 ~]$ **canu \
 -pacbio-raw /usr/local/apps/canu/p6.25x.fastq \
 -usegrid=0 \
 -maxMemory=$(( SLURM\_MEM\_PER\_NODE - 1 )) \
 -maxThreads=$SLURM\_CPUS\_PER\_TASK \
 -p ecoli \
 -d ecoli-auto \
 -genomeSize=4.8m**
  
-- Canu 1.7
--
-- CITATIONS
--
-- Koren S, Walenz BP, Berlin K, Miller JR, Phillippy AM.
-- Canu: scalable and accurate long-read assembly via adaptive k-mer weighting and repeat separation.
-- Genome Res. 2017 May;27(5):722-736.
-- http://doi.org/10.1101/gr.215087.116
[...]


[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. canu.sh). For example:



```

#!/bin/bash
#  this file is called canu_nogrid.sh

cd /data/user/canu
module load canu/1.1

canu \
 -p ecoli -d ecoli-auto \
 -genomeSize=4.8m \
 -pacbio-raw p6.25x.fastq \
    usegrid=0 \
    -maxMemory=$(( SLURM_MEM_PER_NODE - 1 )) \
    -maxThreads=$SLURM_CPUS_PER_TASK \

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] canu.sh
```


Batch jobs using grid options of Canu

In most cases users will want to use the grid options of Canu to distribute the work. Run the Canu command on the
Biowulf login node or in an interactive session, and Canu will submit the jobs appropriately. For example:



```

[user@biowulf canu]$ **module load canu**

[user@biowulf canu]$ **canu -p asm -d lambda -genomeSize=50k -pacbio-raw \
> p6.25x.fastq \
> minReadLength=500 minOverlapLength=500 \
> usegrid=1 \
> gridOptions="--time=30:00 --partition quick" \
> gridOptionsJobName=lam**
-- Detected Java(TM) Runtime Environment '1.8.0_11' (from 'java').
-- Detected 30 CPUs and 118 gigabytes of memory.
-- Detected Slurm with 'sinfo' binary in /usr/local/slurm/bin/sinfo.
--
-- Found  16 hosts with  12 cores and   45 GB memory under Slurm control.
-- Found 594 hosts with  32 cores and  124 GB memory under Slurm control.
-- Found  64 hosts with  32 cores and  124 GB memory under Slurm control.
-- Found  24 hosts with  32 cores and  251 GB memory under Slurm control.
-- Found 384 hosts with  24 cores and   22 GB memory under Slurm control.
-- Found 250 hosts with   8 cores and    6 GB memory under Slurm control.
-- Found 103 hosts with  32 cores and   30 GB memory under Slurm control.
-- Found  16 hosts with  16 cores and   69 GB memory under Slurm control.
-- Found  16 hosts with  16 cores and   69 GB memory under Slurm control.
-- Found 295 hosts with  16 cores and   22 GB memory under Slurm control.
-- Found   4 hosts with  64 cores and 1008 GB memory under Slurm control.
-- Found  64 hosts with  32 cores and   61 GB memory under Slurm control.
-- Found 295 hosts with  32 cores and   61 GB memory under Slurm control.
--
-- Allowed to run under grid control, and use up to   4 compute threads and    3 GB memory for stage 'bogart (unitigger)'.
-- Allowed to run under grid control, and use up to   8 compute threads and    6 GB memory for stage 'mhap (overlapper)'.
-- Allowed to run under grid control, and use up to   8 compute threads and    6 GB memory for stage 'mhap (overlapper)'.
-- Allowed to run under grid control, and use up to   8 compute threads and    6 GB memory for stage 'mhap (overlapper)'.
-- Allowed to run under grid control, and use up to   4 compute threads and    3 GB memory for stage 'read error detection (overlap error adjustment)'.
-- Allowed to run under grid control, and use up to   1 compute thread  and    2 GB memory for stage 'overlap error adjustment'.
-- Allowed to run under grid control, and use up to   4 compute threads and    8 GB memory for stage 'utgcns (consensus'.
-- Allowed to run under grid control, and use up to   1 compute thread  and    2 GB memory for stage 'overlap store parallel bucketizer'.
-- Allowed to run under grid control, and use up to   1 compute thread  and    2 GB memory for stage 'overlapper'.
-- Allowed to run under grid control, and use up to   8 compute threads and    6 GB memory for stage 'overlapper'.
-- Allowed to run under grid control, and use up to   8 compute threads and    6 GB memory for stage 'overlapper'.
-- Allowed to run under grid control, and use up to   4 compute threads and    4 GB memory for stage 'meryl (k-mer counting)'.
-- Allowed to run under grid control, and use up to   4 compute threads and    6 GB memory for stage 'falcon_sense (read correction)'.
-- Allowed to run under grid control, and use up to   8 compute threads and    6 GB memory for stage 'minimap (overlapper)'.
-- Allowed to run under grid control, and use up to   8 compute threads and    6 GB memory for stage 'minimap (overlapper)'.
-- Allowed to run under grid control, and use up to   8 compute threads and    6 GB memory for stage 'minimap (overlapper)'.
----------------------------------------
-- Starting command on Fri Mar 25 14:37:09 2016 with 114.6 GB free disk space

    sbatch \
      --mem=4g \
      --cpus-per-task=1 \
      --time=30:00 \
      --partition quick  \
      -D `pwd` \
      -J "canu_asm_lam" \
      -o /data/user/canu/lambda/canu-scripts/canu.01.out /data/user/canu/lambda/canu-scripts/canu.01.sh
16390160

-- Finished on Fri Mar 25 14:37:09 2016 (lickety-split) with 114.6 GB free disk space
----------------------------------------

```


At various times, the 'sjobs' command will show different Canu jobs running or pending....

```

[user@biowulf canu]$  **sjobs**
User    JobId     JobName       Part   St  Reason      Runtime  Walltime  Nodes  CPUs  Memory    Dependency          Nodelist
================================================================
user  16390181  canu_asm_lam  quick  PD  Dependency     0:00     30:00      1     1  4GB/node  afterany:16390180_*
================================================================

[user@biowulf canu]$ **sjobs**
User    JobId         JobName       Part   St  Reason      Runtime  Walltime  Nodes  CPUs  Memory    Dependency          Nodelist
================================================================
user  16390277_[1]  meryl_asm_la  quick  PD  ---            0:00     30:00      1     4  4GB/node
user  16390278      canu_asm_lam  quick  PD  Dependency     0:00     30:00      1     1  4GB/node  afterany:16390277_*
================================================================

```














