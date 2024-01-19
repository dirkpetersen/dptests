

document.querySelector('title').textContent = 'stampy on Biowulf';
stampy on Biowulf


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



Stampy is a short read mapper for Illumina data with particular emphasis on
mapping reads with a higher number of sequence differences.



### References:


* G. Lunter and M. Goodson. *Stampy: A statistical algorithm for sensitive and fast mapping of Illumina sequence reads*. 
 Genome Reaserch, 2011: 21:936-939.
 [PubMed](https://www.ncbi.nlm.nih.gov/pubmed/20980556)  | 
 [PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3106326/)  | 
 [Journal](http://genome.cshlp.org/content/21/6/936)



Documentation
* [Home page](http://www.well.ox.ac.uk/project-stampy)


Important Notes
* Module Name: stampy (see [the modules page](/apps/modules.html) for more information)
* stampy is a multithreaded application. Please make sure to match the number
 of threads with the number of allocated CPUs.
* Example files in `$STAMPY_TEST_DATA`



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --cpus-per-task=4 --gres=lscratch:10**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144]$ **module load stampy samtools**
[user@cn3144]$ **cd /lscratch/$SLURM\_JOB\_ID**
[user@cn3144]$ **cp $STAMPY\_TEST\_DATA/\* .**
[user@cn3144]$ **ls -lh**
total 68M
-rw-r--r-- 1 user group 56M Jun 21 08:54 ERR458495.fastq.gz
-rwxrwxr-x 1 user group 12M Jun 21 08:54 sacCer3.fa

```

Build a genome and a hash file



```

[user@cn3144]$ **stampy.py -G sacCer3 \
 --species Saccharomyces\_cerevisiae \
 --assembly sacCer3 sacCer3.fa**
stampy: Building genome...
stampy: Input files: ['sacCer3.fa']
stampy: Done
[user@cn3144]$ **stampy.py -g sacCer3 -H sacCer3**
stampy: Building hash table...
stampy: Initializing...
stampy: Counting...
stampy: Initializing hash...         
stampy: Flagging high counts...           
stampy: Creating hash...            
stampy: Writing...                  
stampy: Finished building hash table
stampy: Done

[user@cn3144]$ **ls -lh**
total 87M
-rw-r--r-- 1 user group  56M Jun 21 08:54 ERR458495.fastq.gz
-rwxrwxr-x 1 user group  12M Jun 21 08:54 sacCer3.fa
-rw-r--r-- 1 user group  16M Jun 21 08:58 sacCer3.sthash
-rw-r--r-- 1 user group 2.9M Jun 21 08:58 sacCer3.stidx

```

Align single end data



```

[user@cn3144]$ **stampy.py -t $SLURM\_CPUS\_PER\_TASK -g sacCer3 -h sacCer3 -M ERR458495.fastq.gz \
 | samtools sort -@2 -T ./test -o test.bam**
stampy: Mapping...
stampy: # Nucleotides (all/1/2):        52874862        52874862        0
stampy: # Variants:                     329629  329629  0
stampy: # Fraction:                     0.0062  0.0062  0.0000
stampy: Done

```

End the sinteractive session



```

[user@cn3144]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. stampy.sh) similar to the following example:



```

#! /bin/bash

module load stampy/1.0.31 samtools/1.8 || exit 1
stampy.py -g sacCer3 -h sacCer3 -t $SLURM_CPUS_PER_TASK -M ERR458495.fastq.gz \
  | samtools sort -T /lscratch/$SLURM_JOB_ID/test -@2 -m1g -o test.bam

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=4 --mem=4g --gres=lscratch:10 stampy.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. stampy.swarm). For example:



```

stampy.py -g sacCer3 -h sacCer3 -t $SLURM_CPUS_PER_TASK -M sample1.fastq.gz \
  | samtools sort -T /lscratch/$SLURM_JOB_ID/sample -@2 -m1g -o sample1.bam
stampy.py -g sacCer3 -h sacCer3 -t $SLURM_CPUS_PER_TASK -M sample2.fastq.gz \
  | samtools sort -T /lscratch/$SLURM_JOB_ID/sample -@2 -m1g -o sample2.bam
stampy.py -g sacCer3 -h sacCer3 -t $SLURM_CPUS_PER_TASK -M sample3.fastq.gz \
  | samtools sort -T /lscratch/$SLURM_JOB_ID/sample -@2 -m1g -o sample3.bam

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f stampy.swarm -g 4 -t 4 --gres=lscratch:10 --module stampy/1.0.32,samtools/1.8
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g #  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t #  Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module stampy  Loads the stampy module for each subjob in the swarm 
 | |
 | |
 | |








