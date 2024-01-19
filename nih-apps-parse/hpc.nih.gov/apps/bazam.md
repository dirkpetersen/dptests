

document.querySelector('title').textContent = 'bazam on Biowulf';
bazam on Biowulf


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



A tool to extract paired reads in FASTQ format from coordinate sorted BAM files.




Features
* Bazam is a smarter way to realign reads from one genome to another. If you've tried to use Picard SAMtoFASTQ or samtools bam2fq before and ended up unsatisfied with complicated, long running inefficient pipelines, bazam might be what you wanted. Bazam will output FASTQ in a form that can stream directly into common aligners such as BWA or Bowtie2, so that you can quickly and easily realign reads without extraction to any intermediate format. Bazam can target a specific region of the genome, specified as a region or a gene name if you prefer.




### References:


* Sadedin SP, Oshlack A. *Bazam: a rapid method for read extraction and realignment of high-throughput sequencing data.* Genome Biol. 2019;20(1):78. Published 2019 Apr 18. doi:10.1186/s13059-019-1688-1
PubMed PMID: 30999943; PubMed Central PMCID:PMC6472072.
 [PubMed](https://www.ncbi.nlm.nih.gov/pubmed/30999943) | 
 [PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6472072/) | 
 [Journal](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1688-1)


Documentation
* bazam Main Site:[Main Site](https://github.com/ssadedin/bazam)


Important Notes
* Module Name: bazam (see [the modules page](/apps/modules.html) for more information)
 * Current bazam command lines could be run in two ways:
 
```

	java -jar $BAZAMPATH/bazam.jar --help
	bazam --help
```
* Environment variables set 
	+ $BAZAMPATH* Example files in $BAZAM\_TEST\_DATA



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

[user@cn3144 ~]$ **module load bazam**
[user@cn3144 ~]$ **cp $BAZAM\_TEST\_DATA/\* .**
[user@cn3144 ~]$ **bazam --help**
================================================================================

Bazam

================================================================================

error: Missing required option: bam
usage: java -jar bazam.jar -bam  -L 
 -bam  BAM file to extract read pairs from
 -dr  Specify a read name to debug: processing of the read
 will be verbosey printed
 -f,--filter  Filter using specified groovy expression
 -gene  Extract region of given gene
 -h,--help Show help
 -L,--regions  Regions to include reads (and mates of reads) from
 -n  Concurrency parameter (4)
 -namepos Add original position to the read names
 -o  Output file
 -pad  Amount to pad regions by (0)
 -r1  Output for R1 if extracting FASTQ in separate files
 -r2  Output for R2 if extracting FASTQ in separate files
 -s  Sharding factor: format ,: output only reads
 belonging to shard n of N

This tool is built with Groovy NGS - the Groovy way to work with NGS data.

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Realigning a genome to a new reference using bwa

```

[user@biowulf]$ **sinteractive --cpus-per-task=6 --mem=8g --gres=lscratch:10**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load bazam bwa samtools**
[user@cn3144 ~]$ **cp $BAZAM\_TEST\_DATA/\* .**
[user@cn3144 ~]$ **bazam -bam test.bam | \
bwa mem -t 6 /fdb/igenomes/Homo\_sapiens/UCSC/hg38/Sequence/BWAIndex/genome.fa - | \
samtools view -bhS - | \
samtools sort - -o out.bam -T /lscratch/$SLURM\_JOB\_ID/ -@ 6** 
================================================================================

Bazam

================================================================================

gngs.ToolBase	[1]	INFO	|10:19:19 Auto detected proxy host=dtn05-e0, proxy port=3128
bazam.Bazam	[1]	INFO	|10:19:19 Extracting read pairs from test.bam
bazam.Bazam	[1]	INFO	|10:19:19 Initialising regions to scan from false
gngs.pair.PairScanner	[1]	INFO	|10:19:19 Beginning scan of test.bam
gngs.pair.PairScanner	[1]	INFO	|10:19:19 Created 4 read pair locators
gngs.pair.PairScanner	[1]	INFO	|10:19:20 Stopping parallel threads ...
gngs.pair.PairScanner	[1]	INFO	|10:19:20 Stopping Locator 0
gngs.pair.PairScanner	[1]	INFO	|10:19:20 Stopping Locator 1
gngs.pair.PairScanner	[1]	INFO	|10:19:20 Stopping Locator 2
gngs.pair.PairScanner	[1]	INFO	|10:19:20 Stopping Locator 3
gngs.pair.PairScanner	[1]	INFO	|10:19:20 Stopping Chimeric Locator
gngs.pair.PairScanner	[1]	INFO	|10:19:20 Stopping Formatter
gngs.pair.PairScanner	[1]	INFO	|10:19:20 Stopping Writer
[M::bwa_idx_load_from_disk] read 0 ALT contigs
gngs.pair.PairScanner	[1]	INFO	|10:19:28 Processed 16670 in 9.040 seconds  @ 1843.62/s   (chrX:153649427, loc: 16.7k,16.7k,0 chimeric: 0 formatted: 16.7k, written: 16.7k)
[M::process] read 16670 sequences (1643855 bp)...
[M::mem_process_seqs] Processed 16670 reads in 3.638 CPU sec, 0.612 real sec
[main] Version: 0.7.17-r1188
[main] CMD: bwa mem -t 6 /fdb/igenomes/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/genome.fa -
[main] Real time: 12.340 sec; CPU: 14.019 sec
[bam_sort_core] merging from 0 files and 6 in-memory blocks...
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. bazam.sh). For example:



```

#!/bin/bash
set -e
module load bazam
java -Xmx12g -Dsamjdk.reference_fasta=cram_ref.fasta \
-jar $BAZAMPATH/bazam.jar test.cram >test.fastq.gz
```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=2 --mem=12g bazam.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. bazam.swarm). For example:



```

cd dir1;bazam -bam test1.bam >test1.fastq
cd dir2;bazam -bam test2.bam >test2.fastq
cd dir3;bazam -bam test3.bam >test3.fastq
cd dir4;bazam -bam test4.bam >test4.fastq

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f bazam.swarm [-g #] [-t #] --module bazam
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module bazam Loads the bazam module for each subjob in the swarm
 | |
 | |
 | |








