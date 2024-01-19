

document.querySelector('title').textContent = 'fqtools on Biowulf';
fqtools on Biowulf


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



fqtools are used for manipulating fastq and unaligned bam files. See below for an overview of which tools
are available.



### References:


* A. P. Droop, *fqtools: An efficient software suite for modern FASTQ file manipulation.* 
 Bioinformatics 2016, 32:1883-1884.
 [PubMed](https://www.ncbi.nlm.nih.gov/pubmed/27153699) | 
 [PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4908325/) | 
 [Journal](https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/btw088)


Documentation
* fqtools Main Site: [alastair-droop/fqtools](https://github.com/alastair-droop/fqtools)


Important Notes
* Module Name: fqtools (see [the modules page](/apps/modules.html) for more information)
* Environment variables: $FQTOOLS\_TEST\_DATA



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

[user@cn3144 ~]$ **module load fqtools**
[user@cn3144 ~]$ **fqtools -h**
...
COMMAND:
view      View FASTQ files
head      View the first reads in FASTQ files
count     Count FASTQ file reads
header    View FASTQ file header data
sequence  View FASTQ file sequence data
quality   View FASTQ file quality data
header2   View FASTQ file secondary header data
fasta     Convert FASTQ files to FASTA format
basetab   Tabulate FASTQ base frequencies
qualtab   Tabulate FASTQ quality character frequencies
lengthtab Tabulate FASTQ read lengths
type      Attempt to guess the FASTQ quality encoding type
validate  Validate FASTQ files
find      Find FASTQ reads containing specific sequences
trim      Trim reads in a FASTQ file
qualmap   Translate quality values using a mapping file

[user@cn3144 ~]$ **cp $FQTOOLS\_TEST\_DATA/test.fastq.gz .**

[user@cn3144 ~]$ **fqtools head test.fastq.gz**
@DFXGT8Q1:221:C1F6PACXX:3:1101:1229:2136 1:N:0:
TGTGTTGTCACGCTGCTAATGTCTGCTCTCTCTCGTTTCTTTTTGGAGGC
+
???DDADBC=28AEEEIEDFFI<F@EFIIEEII*??DDEDEIED8@/?;C
...
[user@cn3144 ~]$ **fqtools validate test.fastq.gz**
OK
[user@cn3144 ~]$ **fqtools head -n 20 test.fastq.gz | fqtools quality**
...
[user@cn3144 ~]$ **fqtools type test.fastq.gz**
fastq-sanger
[user@cn3144 ~]$ **fqtools qualtab test.fastq.gz**
!       0
"       0
#       4330007
$       0
%       0
&   59273
'       9373
(       12313
)       79433
...
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. fqtools.sh), which uses the input file 'fqtools.in'. For example:



```

#!/bin/bash
# this file is fqtools.sh
module load fqtools
cp $FQTOOLS_TEST_DATA/test.fastq.gz .
fqtools qualtab test.fastq.gz > test.qualtab
fqtools basetab test.fastq.gz > test.basetab
fqtools lengthtab test.fastq.gz > test.lengthtab
fqtools type test.fastq.gz > test.qual
if [[ "$(fqtools validate test.fastq.gz)" == "OK" ]]; then
    touch test.OK
else
    touch test.FAIL
fi


```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] fqtools.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit set of independent commands requiring identical resources.
Create a swarmfile (e.g. fqtools.swarm). For example:



```

fqtools type sample1.fastq.gz > sample1.qualtype
fqtools type sample2.fastq.gz > sample2.qualtype
fqtools type sample3.fastq.gz > sample3.qualtype

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f fqtools.swarm [-g #] [-t #] --module fqtools
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g #  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t #  Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module fqtools  Loads the fqtools module for each subjob in the swarm 
 | |
 | |
 | |








