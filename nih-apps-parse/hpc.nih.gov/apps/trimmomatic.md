

document.querySelector('title').textContent = 'Trimmomatic on Biowulf';
Trimmomatic on Biowulf


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



Trimmomatic performs a variety of useful trimming tasks for illumina paired-end and single ended data.The selection of trimming steps and their associated parameters are supplied on the command line.

Trimmomatic was developed at the [Usadel lab in Aachen, Germany](http://www.usadellab.org/cms/index.php?page=trimmomatic).

The current trimming steps are:
* ILLUMINACLIP: Cut adapter and other illumina-specific sequences from the read.
* SLIDINGWINDOW: Perform a sliding window trimming, cutting once the average quality within the window falls below a threshold.
* LEADING: Cut bases off the start of a read, if below a threshold quality
* TRAILING: Cut bases off the end of a read, if below a threshold quality
* CROP: Cut the read to a specified length
* HEADCROP: Cut the specified number of bases from the start of the read
* MINLEN: Drop the read if it is below a specified length
* TOPHRED33: Convert quality scores to Phred-33
* TOPHRED64: Convert quality scores to Phred-64



It works with FASTQ (using phred + 33 or phred + 64 quality scores, depending on the Illumina pipeline used), either uncompressed or gzipp'ed FASTQ. Use of gzip format is determined based on the .gz extension.

For single-ended data, one input and one output file are specified, plus the processing steps. For paired-end data, two input files are specified, and 4 output files, 2 for the 'paired' output where both reads survived the processing, and 2 for corresponding 'unpaired' output where a read survived, but the partner read did not.


Use the [modules commands](modules.html) to set up trimmomatic, as
in the example below. The module sets environment variables called 
'TRIMMOJAR' and 'TRIMMOMATIC\_JAR' which point to the location of the
trimmomatic java file. The variable 'TRIMMOMATIC\_JARPATH' points to the 
directory in which the trimmomatic jar file is located.
Fasta files of adapter sequences are included with trimmomatic and can be found at
 **`/usr/local/apps/trimmomatic/&ltversion>/adapters`** 





### References:


* Bolger et al.,Trimmomatic: a flexible trimmer for Illumina sequence data. Bioinformatics (2014). [Link](https://academic.oup.com/bioinformatics/article/30/15/2114/2390096)


Documentation
* [Trimmomatic Main Site](http://www.usadellab.org/cms/index.php?page=trimmomatic)


Important Notes
* Module Name: trimmomatic (see [the modules page](/apps/modules.html) for more information)
* Multithreaded
* environment variables set 
	+ $TRIMMOJAR* Reference data in /usr/local/apps/trimmomatic/&ltversion>/adapters



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem 10g**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load trimmomatic**

[user@cn3144 ~]$ **java -jar $TRIMMOJAR PE -phred33 \**
    **/fdb/app\_testdata/fastq/H\_sapiens/hg100\_1m\_pe1.fq.gz /fdb/app\_testdata/fastq/H\_sapiens/hg100\_1m\_pe2.fq.gz \**
    **output\_forward\_paired.fq.gz output\_forward\_unpaired.fq.gz \**
    **output\_reverse\_paired.fq.gz output\_reverse\_unpaired.fq.gz \**
    **ILLUMINACLIP:/usr/local/apps/trimmomatic/0.39/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36**
TrimmomaticPE: Started with arguments:
 -phred33 /fdb/app_testdata/fastq/H_sapiens/hg100_1m_pe1.fq.gz /fdb/app_testdata/fastq/H_sapiens/hg100_1m_pe2.fq.gz out1 out2 out3 out4 ILLUMINACLIP:/usr/local/apps/trimmomatic/0.39/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
Multiple cores found: Using 2 threads
Using PrefixPair: 'TACACTCTTTCCCTACACGACGCTCTTCCGATCT' and 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT'
ILLUMINACLIP: Using 1 prefix pairs, 0 forward/reverse sequences, 0 forward only sequences, 0 reverse only sequences
Input Read Pairs: 1000000 Both Surviving: 964918 (96.49%) Forward Only Surviving: 25517 (2.55%) Reverse Only Surviving: 7780 (0.78%) Dropped: 1785 (0.18%)
TrimmomaticPE: Completed successfully

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. trimmomatic.sh). For example:



```

#!/bin/bash

ml trimmomatic || exit 1
java -Djava.io.tmpdir=. -jar $TRIMMOJAR PE -phred33 -threads $SLURM_CPUS_PER_TASK \
    SRR292678_1.fastq.gz SRR292678_2.fastq.gz \
    output_forward_paired.fq.gz output_forward_unpaired.fq.gz \
    output_reverse_paired.fq.gz output_reverse_unpaired.fq.gz \
    ILLUMINACLIP:/usr/local/apps/trimmomatic/0.39/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 \
    SLIDINGWINDOW:4:15 MINLEN:36


```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch -c 2 --mem=6g trimmomatic.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. trimmomatic.swarm). For example:



```

java -Djava.io.tmpdir=. -jar $TRIMMOJAR PE -phred33 -threads $SLURM_CPUS_PER_TASK \
    /fdb/app_testdata/fastq/H_sapiens/hg100_1m_pe1.fq.gz /fdb/app_testdata/fastq/H_sapiens/hg100_1m_pe2.fq.gz \
    output_forward_paired.fq.gz output_forward_unpaired.fq.gz \
    output_reverse_paired.fq.gz output_reverse_unpaired.fq.gz \
    ILLUMINACLIP:/usr/local/apps/trimmomatic/0.39/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 \
    SLIDINGWINDOW:4:15 MINLEN:36
[...etc....]

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f trimmomatic.swarm -g 6 -t 8 --module trimmomatic
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module trimmomatic Loads the trimmomatic module for each subjob in the swarm 
 | |
 | |
 | |


















