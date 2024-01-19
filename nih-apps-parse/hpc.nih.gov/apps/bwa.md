

document.querySelector('title').textContent = 'Burrows-Wheeler Alignment (BWA) Tool on Biowulf';
Burrows-Wheeler Alignment (BWA) Tool on Biowulf


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


BWA is a fast light-weighted tool that aligns short sequences to a sequence database, such as the human reference genome. By default, BWA finds an alignment within edit distance 2 to the query sequence, except for disallowing gaps close to the end of the query. It can also be tuned to find a fraction of longer gaps at the cost of speed and of more false alignments.



BWA excels in its speed. Mapping 2 million high-quality 35bp short reads against the human genome can be done in 20 minutes. Usually the speed is gained at the cost of huge memory, disallowing gaps and/or the hard limits on the maximum read length and the maximum mismatches. BWA does not. It is still relatively light-weighted (2.3GB memory for human alignment), performs gapped alignment, and does not set a hard limit on read length or maximum mismatches.



Given a database file in FASTA format, BWA first builds BWT index with the 'index' command. The alignments in suffix array (SA) coordinates are then generated with the 'aln' command. The resulting file contains ALL the alignments found by BWA. The 'samse/sampe' command converts SA coordinates to chromosomal coordinates. For single-end reads, most of computing time is spent on finding the SA coordinates (the aln command). For paired-end reads, half of computing time may be spent on pairing (the sampe command) given 32bp reads. Using longer reads would reduce the fraction of time spent on pairing because each end in a pair would be mapped to fewer places.



### References:


* [Li H. and Durbin R. (2009) Fast and accurate short read alignment with Burrows-Wheeler Transform. Bioinformatics, 25:1754-60.](http://www.ncbi.nlm.nih.gov/pubmed/19451168)* If you use BWA-SW, please cite:  

[Li H. and Durbin R. (2010) Fast and accurate long-read alignment with Burrows-Wheeler Transform. Bioinformatics, Epub.](http://www.ncbi.nlm.nih.gov/pubmed/20080505)


Documentation
* [bwa Main Site](http://bio-bwa.sourceforge.net/)
* Type **bwa** at the prompt
* Type **bwa *command*** at the prompt, e.g. **bwa index**


Important Notes
* Module Name: bwa (see [the modules page](/apps/modules.html) for more information)
* Multithreaded
* Example files in **/fdb/app\_testdata/fastq/**
* Reference data in:
	+ **/fdb/bwa/**
	+ **/fdb/igenomes/*organism*/*source*/*build*/Sequence/BWAIndex/**, where  
	
		- ***organism*** is the specific organism of interest (Gallus\_gallus, Rattus\_norvegicus, etc.)
		- ***source*** is the source for the sequence (NCBI, Ensembl, UCSC)
		- ***build*** is the specific genome draft of interest (hg19, build37.2, GRCh37)



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --cpus-per-task=8 --mem=24g**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ ml bwa
[user@cn3144 ~]$ bwa aln -t 8 /fdb/bwa/indexes/hg19.fa file.fastq > file.sai

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. bwa.sh). For example:



```

module load bwa
cd /data/$USER/dir
bwa index -a bwtsw file.csfasta
bwa aln -t $SLURM_CPUS_PER_TASK file.csfasta file.fastq > file.sai
bwa samse file.csfasta file.sai file.fastq > file.sam

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] bwa.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. bwa.swarm). For example:



```

cd /data/$USER/dir1; bwa aln -t $SLURM_CPUS_PER_TASK file.csfasta file.fastq > file.sai
cd /data/$USER/dir2; bwa aln -t $SLURM_CPUS_PER_TASK file.csfasta file.fastq > file.sai
cd /data/$USER/dir3; bwa aln -t $SLURM_CPUS_PER_TASK file.csfasta file.fastq > file.sai

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f bwa.swarm [-g #] [-t #] --module bwa
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module bwa Loads the bwa module for each subjob in the swarm 
 | |
 | |
 | |








