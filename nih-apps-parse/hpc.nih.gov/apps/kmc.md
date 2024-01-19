

document.querySelector('title').textContent = 'kmc on Biowulf';
kmc on Biowulf


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



KMC is a program to create and access databases for counting k-mers from
fastq or fasta files.



Documentation
* [Home page](http://sun.aei.polsl.pl/kmc)
* [GitHub](https://github.com/marekkokot/KMC)


Important Notes
* Module Name: kmc (see [the modules page](/apps/modules.html) for more information)
* kmc is a multithreaded application. Make sure to match allocated CPUs with the number of threads
* Example files in `$KMC_TEST_DATA`



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --gres=lscratch:10 --mem=10g --cpus-per-task=2**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144]$ **module load kmc**
[user@cn3144]$ **kmc**
K-Mer Counter (KMC) ver. 3.0.0 (2017-01-28)
Usage:
 kmc [options] <input_file_name> <output_file_name> <working_directory>
 kmc [options] <@input_file_names> <output_file_name> <working_directory>
Parameters:
  input_file_name - single file in FASTQ format (gziped or not)
  @input_file_names - file name with list of input files in FASTQ format (gziped or not)
Options:
  -v - verbose mode (shows all parameter settings); default: false
  -k<len> - k-mer length (k from 1 to 256; default: 25)
  -m<size> - max amount of RAM in GB (from 1 to 1024); default: 12
  -sm - use strict memory mode (memory limit from -m<n> switch will not be exceeded)
  -p<par> - signature length (5, 6, 7, 8, 9, 10, 11); default: 9
  -f<a/q/m> - input in FASTA format (-fa), FASTQ format (-fq) or multi FASTA (-fm); default: FASTQ
  -ci<value> - exclude k-mers occurring less than <value> times (default: 2)
  -cs<value> - maximal value of a counter (default: 255)
  -cx<value> - exclude k-mers occurring more of than <value> times (default: 1e9)
  -b - turn off transformation of k-mers into canonical form
  -r - turn on RAM-only mode
  -n<value> - number of bins
  -t<value> - total number of threads (default: no. of CPU cores)
  -sf<value> - number of FASTQ reading threads
  -sp<value> - number of splitting threads
  -sr<value> - number of threads for 2nd stage
Example:
kmc -k27 -m24 NA19238.fastq NA.res \data\kmc_tmp_dir\
kmc -k27 -m24 @files.lst NA.res \data\kmc_tmp_dir\

[user@cn3144]$ **cd /lscratch/$SLURM\_JOB\_ID**
[user@cn3144]$ **cp $KMC\_TEST\_DATA/ENCFF001KPB.fastq.gz .**
[user@cn3144]$ **mkdir ENCFF001KPB.tmp**
[user@cn3144]$ **kmc -t2 ENCFF001KPB.fastq.gz ENCFF001KPB.kmc ENCFF001KPB.tmp**
*****************
Stage 1: 100%
Stage 2: 100%
1st stage: 8.87882s
2nd stage: 6.6567s
Total    : 15.5355s
Tmp size : 173MB

Stats:
   No. of k-mers below min. threshold :     66106376
   No. of k-mers above max. threshold :            0
   No. of unique k-mers               :     79254698
   No. of unique counted k-mers       :     13148322
   Total no. of k-mers                :    108850013
   Total no. of reads                 :      9157799
   Total no. of super-k-mers          :     19752822


[user@cn3144]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. kmc.sh) similar to the following example:



```

#! /bin/bash
set -e

tmp=/lscratch/${SLURM_JOB_ID}
module load kmc/3.0.0 || exit 1

kmc -t$(( SLURM_CPUS_PER_TASK - 2 )) -m15 -sm \
  $KMC_TEST_DATA/ENCFF001KPB.fastq.gz \
  ENCFF001KPB.kmer \
  ${tmp}

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=4 --mem=10g --gres=lscratch:10 kmc.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. kmc.swarm). For example:



```

tmp=/lscratch/${SLURM_JOB_ID} \
  && kmc -t$(( SLURM_CPUS_PER_TASK - 2 )) -m15 -sm 1.fastq.gz 1.kmer ${tmp}
tmp=/lscratch/${SLURM_JOB_ID} \
  && kmc -t$(( SLURM_CPUS_PER_TASK - 2 )) -m15 -sm 2.fastq.gz 2.kmer ${tmp}

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f kmc.swarm [-g #] [-t #] --module kmc/3.0.0 --gres=lscratch:10
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g #  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t #  Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module kmc  Loads the kmc module for each subjob in the swarm 
 | |
 | |
 | |








