

document.querySelector('title').textContent = 'winnowmap on Biowulf';
winnowmap on Biowulf


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



From the documentation:



>  Winnowmap is a long-read mapping algorithm optimized for mapping ONT and
>  PacBio reads to repetitive reference sequences. Winnowmap development began on
>  top of minimap2 codebase, and since then we have incorporated the following two
>  ideas to improve mapping accuracy within repeats. ... 


### References:


* Chirag Jain, Arang Rhie, Nancy Hansen, Sergey Koren and Adam Phillippy. 
 *A long read mapping method for highly repetitive reference sequences*. 
 [BioRxiv, 2020.](https://www.biorxiv.org/content/10.1101/2020.11.01.363887v1)



Documentation
* winnowmap on [GitHub](https://github.com/marbl/Winnowmap)


Important Notes
* Module Name: winnowmap (see [the modules page](/apps/modules.html) for more information)
* winnowmap is a multithreaded application. Make sure to match the number of threads with the number of allocated CPUs
* Example files in `$WINNOWMAP_TEST_DATA`



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --cpus-per-task=6 --mem=8g --gres=lscratch:10**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144]$ **cd /lscratch/$SLURM\_JOB\_ID**
[user@cn3144]$ **module load winnowmap**
[user@cn3144]$ **cp -rL ${WINNOWMAP\_TEST\_DATA:-none}/\* .**
[user@cn3144]$ ## note that memory is slighly less the amount of allocated memory
[user@cn3144]$ **meryl count k=15 threads=$SLURM\_CPUS\_PER\_TASK memory=7 \
 output merylDB Zymo-Isolates-SPAdes-Illumina.fasta**

Found 1 command tree.

Counting 68 (estimated) million canonical 15-mers from 1 input file:
    sequence-file: Zymo-Isolates-SPAdes-Illumina.fasta


SIMPLE MODE
-----------

  15-mers
    -> 1073741824 entries for counts up to 65535.
    -> 16 Gbits memory used

  71433493 input bases
    -> expected max count of 285733, needing 4 extra bits.
    -> 4096 Mbits memory used

  2560 MB memory needed


COMPLEX MODE
------------

prefix     # of   struct   kmers/    segs/      min     data    total
  bits   prefix   memory   prefix   prefix   memory   memory   memory
------  -------  -------  -------  -------  -------  -------  -------
     1     2  P    37 kB    34 MM  1976  S   128 kB   247 MB   247 MB
     2     4  P    42 kB    17 MM   954  S   256 kB   238 MB   238 MB
     3     8  P    54 kB  8719 kM   460  S   512 kB   230 MB   230 MB
     4    16  P    78 kB  4359 kM   222  S  1024 kB   222 MB   222 MB
     5    32  P   127 kB  2179 kM   107  S  2048 kB   214 MB   214 MB
     6    64  P   228 kB  1089 kM    52  S  4096 kB   208 MB   208 MB
     7   128  P   429 kB   544 kM    25  S  8192 kB   200 MB   200 MB
     8   256  P   832 kB   272 kM    12  S    16 MB   192 MB   192 MB
     9   512  P  1640 kB   136 kM     6  S    32 MB   192 MB   193 MB
    10  1024  P  3256 kB    68 kM     3  S    64 MB   192 MB   195 MB  Best Value!
    11  2048  P  6496 kB    34 kM     2  S   128 MB   256 MB   262 MB
    12  4096  P    12 MB    17 kM     1  S   256 MB   256 MB   268 MB
    13  8192  P    25 MB  8720  M     1  S   512 MB   512 MB   537 MB
    14    16 kP    50 MB  4360  M     1  S  1024 MB  1024 MB  1074 MB
    15    32 kP   101 MB  2180  M     1  S  2048 MB  2048 MB  2149 MB
    16    64 kP   202 MB  1090  M     1  S  4096 MB  4096 MB  4298 MB


FINAL CONFIGURATION
-------------------

Configured complex mode for 0.191 GB memory per batch, and up to 1 batch.

Start counting with THREADED method.
Used 0.285 GB out of 125.715 GB to store      1048320 kmers.
Used 0.410 GB out of 125.715 GB to store     54504240 kmers.

Writing results to 'merylDB', using 4 threads.
finishIteration()--

Finished counting.

Cleaning up.

Bye.
[user@cn3144]$ **meryl print greater-than distinct=0.9998 merylDB > repetitive\_k15.txt**

Found 1 command tree.

PROCESSING TREE #1 using 1 thread.
  opGreaterThan
    threshold=12
    fraction-distinct=0.999800
    merylDB
    print to (stdout)

Cleaning up.

Bye.
[user@cn3144]$ **ls -lh**
total 715M
drwxr-xr-x 2 user group  12K Apr  2 08:30 merylDB
-rw-r--r-- 1 user group 647M Apr  2 08:27 ont_zymo.fastq.gz
-rw-r--r-- 1 user group 198K Apr  2 08:31 repetitive_k15.txt
-rw-r--r-- 1 user group  69M Apr  2 08:27 Zymo-Isolates-SPAdes-Illumina.fasta
[user@cn3144]$ **head repetitive\_k15.txt**
AAAAAAAAAAAAAAA 8426
AAAAAAAAAAAAAAC 241
AAAAAAAAAAAAAAT 307
AAAAAAAAAAAAAAG 365
AAAAAAAAAAAAACA 143
AAAAAAAAAAAAACC 69
AAAAAAAAAAAAACT 72
AAAAAAAAAAAAACG 44
AAAAAAAAAAAAATA 154
AAAAAAAAAAAAATC 44
[user@cn3144]$ **module load samtools**
[user@cn3144]$ **winnowmap -W repetitive\_k15.txt -t $SLURM\_CPUS\_PER\_TASK -ax map-ont \
 Zymo-Isolates-SPAdes-Illumina.fasta ont\_zymo.fastq.gz \
 | samtools sort -T /lscratch/$SLURM\_JOB\_ID/part -@2 -m1G -O BAM -o ont\_zymo.bam**
[M::mm_idx_gen::0.000*32.22] reading downweighted kmers
[M::mm_idx_gen::0.004*1.69] collected downweighted kmers, no. of kmers read=10651
[M::mm_idx_gen::0.004*1.69] saved the kmers in a bloom filter: hash functions=2 and size=153136
[M::mm_idx_gen::4.113*1.06] collected minimizers
[M::mm_idx_gen::4.251*1.14] sorted minimizers
[M::main::4.251*1.14] loaded/built the index for 2428 target sequence(s)
[M::main::4.251*1.14] running winnowmap in SV-aware mode
[M::main::4.251*1.14] stage1-specific parameters minP:1000, incP:4.00, maxP:16000, sample:1000, min-qlen:10000, min-qcov:0.5, min-mapq:5, mid-occ:5000
[M::main::4.251*1.14] stage2-specific parameters s2_bw:2000, s2_zdropinv:25
[M::mm_idx_stat] kmer size: 15; skip: 50; is_hpc: 0; #seq: 2428
[M::mm_idx_stat::4.287*1.14] distinct minimizers: 2443036 (89.73% are singletons); average occurrences: 1.125; average spacing: 25.530
[M::worker_pipeline::135.789*4.94] mapped 160000 sequences
[M::main] Version: 2.0, pthreads=6, omp_threads=3
[M::main] CMD: winnowmap -W repetitive_k15.txt -t 6 -ax map-ont Zymo-Isolates-SPAdes-Illumina.fasta ont_zymo.fastq.gz
[M::main] Real time: 135.845 sec; CPU: 670.800 sec; Peak RSS: 3.423 GB

[user@cn3144]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. winnowmap.sh). For example, for aligning reads
to a reference with known repetitive k-mers:



```

#!/bin/bash
module load winnowmap/2.0 samtools/1.11
winnowmap -W repetitive_k15.txt -t $SLURM_CPUS_PER_TASK -ax map-ont \
    ${WINNOWMAP_TEST_DATA:-none}/Zymo-Isolates-SPAdes-Illumina.fasta \
    ${WINNOWMAP_TEST_DATA:-none}/ont_zymo.fastq.gz \
| samtools sort -T /lscratch/$SLURM_JOB_ID/part -@2 -m1G -O BAM -o ont_zymo.bam

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=6 --mem=8g --gres=lscratch:10 winnowmap.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. winnowmap.swarm). For example:



```

winnowmap -W repetitive_k15.txt -t $SLURM_CPUS_PER_TASK -ax map-ont \
    genome.fasta sample1.fastq.gz \
| samtools sort -T /lscratch/$SLURM_JOB_ID/part -@2 -m1G -O BAM sample1.bam
winnowmap -W repetitive_k15.txt -t $SLURM_CPUS_PER_TASK -ax map-ont \
    genome.fasta sample2.fastq.gz \
| samtools sort -T /lscratch/$SLURM_JOB_ID/part -@2 -m1G -O BAM sample2.bam

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f winnowmap.swarm -g 8 -t 6 --module winnowmap/2.0,samtools/1.11
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g #  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t #  Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module winnowmap  Loads the winnowmap module for each subjob in the swarm 
 | |
 | |
 | |








