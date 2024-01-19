

document.querySelector('title').textContent = 'fastp on Biowulf';
fastp on Biowulf


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



A tool designed to provide fast all-in-one preprocessing for FastQ files.
This tool is developed in C++ with multithreading supported to afford high performance.




Features
* 0.comprehensive quality profiling for both before and after filtering data (quality curves, base contents, KMER, Q20/Q30, GC Ratio, duplication, adapter contents...)
 * 1.filter out bad reads (too low quality, too short, or too many N...)
 * 2.cut low quality bases for per read in its 5' and 3' by evaluating the mean quality from a sliding window (like Trimmomatic but faster).
 * 3.trim all reads in front and tail
 * 4.cut adapters. Adapter sequences can be automatically detected, which means you don't have to input the adapter sequences to trim them.
 * 5.correct mismatched base pairs in overlapped regions of paired end reads, if one base is with high quality while the other is with ultra low quality
 * 6.trim polyG in 3' ends, which is commonly seen in NovaSeq/NextSeq data. Trim polyX in 3' ends to remove unwanted polyX tailing (i.e. polyA tailing for mRNA-Seq data)
 * 7.preprocess unique molecular identifier (UMI) enabled data, shift UMI to sequence name.
 * 8.report JSON format result for further interpreting.
 * 9.visualize quality control and filtering results on a single HTML page (like FASTQC but faster and more informative).
 * 10.split the output to multiple files (0001.R1.gz, 0002.R1.gz...) to support parallel processing. Two modes can be used, limiting the total split file number, or limitting the lines of each split file.
 * 11.support long reads (data from PacBio / Nanopore devices).
 * 12.support reading from STDIN and writing to STDOUT
 * 13.support interleaved input




### References:


* Chen S, Zhou Y, Chen Y, Gu J. *fastp: an ultra-fast all-in-one FASTQ
preprocessor*. Bioinformatics. 2018 Sep 1;34(17):i884-i890. doi:
10.1093/bioinformatics/bty560. PubMed PMID: 30423086; PubMed Central PMCID:
PMC6129281.
 [PubMed](https://www.ncbi.nlm.nih.gov/pubmed/30423086) | 
 [PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6129281/) | 
 [Journal](https://academic.oup.com/bioinformatics/article/34/17/i884/5093234)


Documentation
* fastp Main Site:[Main Site](https://github.com/OpenGene/fastp)


Important Notes
* Module Name: fastp (see [the modules page](/apps/modules.html) for more information)
 * Multithreaded
 * This application produces HTML reports. You can use [hpcdrive](/docs/hpcdrive.html) to view these reports on your local workstation.
* Environment variables set 
	+ $FASTP\_TEST\_DATA* Example files in $FASTP\_TEST\_DATA



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

[user@cn3144 ~]$ **module load fastp**
[user@cn3144 ~]$ **cp $FASTP\_TEST\_DATA/\* .**
[user@cn3144 ~]$ **fastp --in1 R1.fq --in2 R2.fq**
Read1 before filtering:
total reads: 9
total bases: 1208
Q20 bases: 1078(89.2384%)
Q30 bases: 1005(83.1954%)

Read2 before filtering:
total reads: 9
total bases: 1359
Q20 bases: 1100(80.9419%)
Q30 bases: 959(70.5666%)

Read1 after filtering:
total reads: 8
total bases: 1208
Q20 bases: 1078(89.2384%)
Q30 bases: 1005(83.1954%)

Read2 aftering filtering:
total reads: 8
total bases: 1208
Q20 bases: 991(82.0364%)
Q30 bases: 874(72.351%)

Filtering result:
reads passed filter: 16
reads failed due to low quality: 0
reads failed due to too many N: 0
reads failed due to too short: 2
reads with adapter trimmed: 0
bases trimmed due to adapters: 0

Duplication rate: 62.5%

Insert size peak (evaluated by paired-end reads): 187

JSON report: fastp.json
HTML report: fastp.html

fastp --in1 R1.fq --in2 R2.fq
fastp v0.20.1, time used: 0 seconds
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. fastp.sh). For example:



```

#!/bin/bash
set -e
module load fastp
fastp --in1 R1.fq --in2 R2.fq

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=12 --mem=2g fastp.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. fastp.swarm). For example:



```

cd dir1;fastp --in1 R1.fq --html R1.html
cd dir2;fastp --in1 R2.fq --html R2.html
cd dir3;fastp --in1 R3.fq --html R3.html
cd dir4;fastp --in1 R4.fq --html R4.html

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f fastp.swarm [-g #] [-t #] --module fastp
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module fastp Loads the fastp module for each subjob in the swarm
 | |
 | |
 | |








