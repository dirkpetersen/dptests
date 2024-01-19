

document.querySelector('title').textContent = 'bismark on Biowulf';
bismark on Biowulf


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



Bismark is a program to map bisulfite treated sequencing reads to a genome of interest and perform methylation calls in a single step. The output can be easily imported into a genome viewer, such as SeqMonk, and enables a researcher to analyse the methylation levels of their samples straight away. It's main features are:

* Bisulfite mapping and methylation calling in one single step
* Supports single-end and paired-end read alignments
* Supports ungapped and gapped alignments
* Alignment seed length, number of mismatches etc. are adjustable
* Output discriminates between cytosine methylation in CpG, CHG and CHH context





### References:


* [Krueger, Felix, and Simon R. Andrews. "Bismark: a flexible aligner and methylation caller for Bisulfite-Seq applications." *Bioinformatics* 27.11 (2011): 1571-1572.](http://bioinformatics.oxfordjournals.org/content/27/11/1571.short)


Documentation
Bismark is extensively documented. To read the help doc, type bismark --help. See also:
* [Home page at Babraham Bioinformatics](http://www.bioinformatics.babraham.ac.uk/projects/bismark/)
* [GitHub repository](https://github.com/FelixKrueger/Bismark)
* [Bismark User Guide](http://www.bioinformatics.babraham.ac.uk/projects/bismark/Bismark_User_Guide.pdf)


Important Notes
* Module Name: bismark (see [the modules page](/apps/modules.html) for more information)
* Multithreaded
* environment variables set 
	+ BISMARK\_HOME* Example file: $BISMARK\_HOME/test\_data.fastq* Reference data in /fdb/bismark/



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
 The first step in using Bismark is to prepare a genome using the bismark\_genome\_preparation command. This creates a subdirectory in the fasta file directory called Bisulfite\_Genome. To prevent unnecessary data duplication, the human genome (hg19 and hg38) has been processed for use with the bowtie2 aligner. It resides under /fdb/bismark. If you would like other genomes to be prepared please contact staff@hpc.nih.gov.

A fastq file with test data is provided by the authors. In this example a user copies this testfile along with X and Y chromosome data from the human genome to thier local space and runs though the basic analysis steps. (User input in bold)

Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=6g**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job


[user@cn3144 ~]$ **bismark --help**

     This program is free software: you can redistribute it and/or modify
     it under the terms of the GNU General Public License as published by
     the Free Software Foundation, either version 3 of the License, or
     (at your option) any later version.

     This program is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
...

```

[Full output from the help command](/apps/bismark/bismark_help.html)    



```

[user@cn3144 ~]$ **mkdir -p /data/$USER/bismark\_test/XandY**
[user@cn3144 ~]$ **cd /data/$USER/bismark\_test**
[user@cn3144 bismark_test]$ **cp /fdb/genome/human-feb2009/chrX.fa ./XandY**
[user@cn3144 bismark_test]$ **cp /fdb/genome/human-feb2009/chrY.fa ./XandY**
[user@cn3144 bismark_test]$ **cp $BISMARK\_HOME/test\_data.fastq .**
[user@cn3144 bismark_test]$ **bismark\_genome\_preparation XandY**
Writing bisulfite genomes out into a single MFA (multi FastA) file
Bisulfite Genome Indexer version v0.16.0 (last modified 25 August 2015)
Step I - Prepare genome folders - completed
Total number of conversions performed:
[....]

```

[Full output from the prep command](/apps/bismark/bismark_prep_output.html)    



```

[user@cn3144 ~]$ **bismark XandY test\_data.fastq**
Path to Bowtie 2 specified as: bowtie2
Output format is BAM (default)
Alignments will be written out in BAM format. Samtools found here: '/usr/local/apps/samtools/1.3.1/bin/samtools'
Reference genome folder provided is XandY/	(absolute path is '/spin1/users/user/bismark_test/XandY/)'
FastQ format assumed (by default)
Files to be analysed:
test_data.fastq
Library is assumed to be strand-specific (directional), alignments to strands complementary to the original top or bottom strands will be ignored (i.e. not performed!)
...

```

[Full output from the run command](/apps/bismark/bismark_run_output.html)    



```

[user@cn3144 bismark_test]$ **bismark\_methylation\_extractor test\_data\_bismark\_bt2.bam**

 *** Bismark methylation extractor version v0.16.0 ***

Trying to determine the type of mapping from the SAM header line of file test_data_bismark_bt2.bam
Treating file(s) as single-end data (as extracted from @PG line)

Setting core usage to single-threaded (default). Consider using --multicore &ltint> to speed up the extraction process.

Summarising Bismark methylation extractor parameters:
===============================================================

...

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


[Full output from the extract command](/apps/bismark/bismark_extract_output.html)    


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. bismark.sh). For example:



```

#!/bin/bash
# this file is called bismark.sh
set -e

module load bismark
cd /data/$USER/bismark_test
bismark_genome_preparation XandY
bismark XandY test_data.fastq
bismark_methylation_extractor test_data_bismark_bt2.bam

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch bismark.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.

One strategy for running a swarm of bismark jobs would be to set up multiple directories containing gene sequences. Once you have done that, you can set up a swarm command file containing one line for each of your bismark runs. 
Sample swarm command file



```

# --------file myjobs.swarm----------
bismark directory1 test_data.fastq
bismark directory2 test_data.fastq
bismark directory3 test_data.fastq
....
bismark directoryN test_data.fastq
# -----------------------------------

```

Submit this set of runs to the batch system by typing



```

[user@biowulf ~]$ **swarm --module bismark -f myjobs.swarm**

```

















