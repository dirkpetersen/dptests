

document.querySelector('title').textContent = 'cutadapt on Biowulf';
cutadapt on Biowulf


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


Cutadapt removes adapter sequences, primers, poly-A tails, low quality
segments, and other unwanted sequence from your high-throughput sequencing 
reads.


### References:


* Marcel Martin. *Cutadapt removes adapter sequences from high-
 throughput sequencing reads.*[DOI:10.14806/ej.17.1.200](http://dx.doi.org/10.14806/ej.17.1.200)


Documentation
* [GitHub](https://github.com/marcelm/cutadapt)
* [Documentation](https://cutadapt.readthedocs.org/en/stable/)


Important Notes
* Module Name: cutadapt (see [the modules page](/apps/modules.html) for more information)
* Starting with version 1.15 cutadapt can make use or multiple cores. Note that in future 
 versions cutadapt may attempt to autodetect cores which may lead to overloading Slurm jobs. That means the
 number of cores should be set explicitly.
* Example data can be found in `$CUTADAPT_TEST_DATA`



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

[user@cn3144 ~]$ **module load cutadapt**
[user@cn3144 ~]$ # copy some paired end RNASeq data (Illumina)
[user@cn3144 ~]$ **cp $CUTADAPT\_TEST\_DATA/\* .**
[user@cn3144 ~]$ **ls -lh**
-rw-r--r-- 1 user group 45M Feb  2 07:22 read1_1000k.fastq.gz
-rw-r--r-- 1 user group 35M Feb  2 07:22 read2_1000k.fastq.gz

[user@cn3144 ~]$ **cutadapt -q 10 --minimum-length 25 \
 -a AGATCGGAAGAGC -A AGATCGGAAGAGC \
 -o read1\_trimmed.fastq.gz -p read2\_trimmed.fastq.gz \
 read1\_1000k.fastq.gz read2\_1000k.fastq.gz**
This is cutadapt 1.15 with Python 3.6.4
Command line parameters: -q 10 --minimum-length 25 -a AGATCGGAAGAGC -A AGATCGGAAGAGC -o read1.fastq.gz -p read2.fastq gz read1_1000k.fastq.gz read2_1000k.fastq.gz
Running on 1 core
Trimming 2 adapters with at most 10.0% errors in paired-end mode ...
Finished in 53.95 s (54 us/read; 1.11 M reads/minute).

=== Summary ===

Total read pairs processed:          1,000,000
  Read 1 with adapter:                  22,329 (2.2%)
  Read 2 with adapter:                  23,325 (2.3%)
Pairs that were too short:             221,834 (22.2%)
Pairs written (passing filters):       778,166 (77.8%)

Total basepairs processed:   100,000,000 bp
  Read 1:    50,000,000 bp
  Read 2:    50,000,000 bp
Quality-trimmed:              29,853,878 bp (29.9%)
  Read 1:    14,926,939 bp
  Read 2:    14,926,939 bp
Total written (filtered):     77,202,392 bp (77.2%)
  Read 1:    38,478,257 bp
  Read 2:    38,724,135 bp

=== First read: Adapter 1 ===

Sequence: AGATCGGAAGAGC; Type: regular 3'; Length: 13; Trimmed: 22329 times.

No. of allowed errors:
0-9 bp: 0; 10-13 bp: 1

Bases preceding removed adapters:
  A: 26.6%
  C: 28.3%
  G: 28.1%
  T: 16.9%
  none/other: 0.1%

Overview of removed sequences
length  count   expect  max.err error counts
3       17956   15625.0 0       17956
4       3301    3906.2  0       3301
...

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. cutadapt.sh), which uses the input file 'cutadapt.in'. For example:



```

#! /bin/bash
set -e

r1=fastq/read1.fastq.gz
r2=fastq/read2.fastq.gz

module load cutadapt/1.15 || exit 1
cutadapt -q 10 --trim-n --minimum-length 25 \
  --cores=$SLURM_CPUS_PER_TASK \
  -a AGATCGGAAGAGC -A AGATCGGAAGAGC \
  -o fastq_clean/${r1#fastq} -p fastq_clean/${r2#fastq} \
  $r1 $r2

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=4 [--mem=#] cutadapt.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. cutadapt.swarm). For example:



```

cutadapt -q10, --trim-n --minimum-length 25 -a AGATCGGAAGAGC -o clean1.fq.gz dirty1.fq.gz
cutadapt -q10, --trim-n --minimum-length 25 -a AGATCGGAAGAGC -o clean2.fq.gz dirty2.fq.gz
cutadapt -q10, --trim-n --minimum-length 25 -a AGATCGGAAGAGC -o clean2.fq.gz dirty2.fq.gz

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f cutadapt.swarm [-g #] [-t #] --module cutadapt
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g #  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t #  Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module cutadapt  Loads the cutadapt module for each subjob in the swarm 
 | |
 | |
 | |








