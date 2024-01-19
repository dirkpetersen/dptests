

document.querySelector('title').textContent = 'tophat on Biowulf';
tophat on Biowulf


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


TopHat is a splice junction mapper for RNA-Seq reads. It aligns RNA-Seq reads
to mammalian-sized genomes using the short read aligner Bowtie, and then
analyzes the mapping results to identify splice junctions between exons.


  

  

### References:


* Cole Trapnell, L. Pachter L, and S. L. Salzberg.
 [**TopHat: discovering splice junctions with RNA-Seq.**](http://www.ncbi.nlm.nih.gov/pubmed/19289445)
*Bioinformatics 2009, 25:1105-1111** Daehwan Kim, G. Pertea, C. Trapnell, H. Pimentel, R. Kelley, S. L. Salzberg SL. 
 [**TopHat2: accurate alignment of transcriptomes in the presence of insertions, deletions and gene fusions.**](http://www.ncbi.nlm.nih.gov/pubmed/23618408)
*Genome Biology 2013, 14:R36.** Daehwan Kim, and S. L. Salzberg.
 [**TopHat-Fusion: an algorithm for discovery of novel fusion transcripts.**](http://www.ncbi.nlm.nih.gov/pubmed/21835007)
*Genome Biology 2011, 12:R72.*



Please note that TopHat has entered a low maintenance, low support stage as it is now largely superseded by [HISAT2](hisat.html) which provides the same core functionality (i.e. spliced alignment of RNA-Seq reads), in a more accurate and much more efficient way.



Documentation
* [Home page](http://ccb.jhu.edu/software/tophat/index.shtml)
* [Manual](http://ccb.jhu.edu/software/tophat/manual.shtml)
* [Google group](https://groups.google.com/forum/#!forum/tuxedo-tools-users)


Important Notes
* Module Name: tophat (see [the modules page](/apps/modules.html) for more information)
* Multithreaded
* Environment variables set 
	+ PATH


TopHat makes use of either `bowtie2` (the default) or `bowtie` (optional; 
necessary for colorspace reads as bowtie2 does not support colorspace). TopHat
is limited to a maximal read length of 1024nts.


TopHat2 now includes TopHat-Fusion's ability to look for fusions between
different transcripts (`--fusion-search`).



Index files


TopHat will need either `bowtie2` or `bowtie` indices which are available
as part of the igenomes package at


**`/fdb/igenomes/[organism]/[source]/[build]/Sequence/Bowtie[2]Index/*`**


* **`[organism]`** is the specific organism of interest
 (Gallus\_gallus, Rattus\_norvegicus, etc.)
* **`[source]`** is the source for the sequence (NCBI,
 Ensembl, UCSC)
* **`[build]`** is the specific genome build of interest
 (hg19, build37.2, GRCh37)


More information on the locally available igenomes builds/organisms
is available from our [scientific database index](https://hpc.cit.nih.gov/apps/db.php?f=Igenomes).
For more information about igenomes in general,
[iGenomes readme](ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/README.txt).



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive --cpus-per-task=8 --mem=10g**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **cd /data/$USER/test\_data**
[user@cn3144 ~]$ **module load tophat**
[user@cn3144 ~]$ **tophat -o job1 -p8 /path/to/genome /path/to/fastq.gz**
...
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. tophat.sh). For example:



```

#! /bin/bash
module load tophat
cd /data/${USER}/test_data
tophat \
  --output-dir=./tophat_test \
  --min-anchor-length=10 \
  --num-threads=$SLURM_CPUS_PER_TASK \
  --b2-sensitive \
  /fdb/igenomes/Mus_musculus/UCSC/mm9/Sequence/Bowtie2Index/genome \
  fastq/rnaseq_500k.fastq.gz

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=16 --mem=20g tophat.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Note that tophat by default creates
a ouput directory called `tophat_out`. Therefore, to run multiple jobs in
parallel, the jobs either have to be run in separate directories or
the command line has to specify the `-o / --output_dir` for each job.


Create a swarmfile (e.g. tophat.swarm). For example:



```

tophat -o job1 -p ${SLURM_CPUS_PER_TASK} --b2-sensitive \
  --GTF=annot/140609_refseq_nomir_nosnor.gtf --no-novel-juncs \    
  /fdb/igenomes/Mus_musculus/UCSC/mm9/Sequence/Bowtie2Index/genome \
  fastq/rnaseq_500k.fastq.gz
tophat -o job2 -p ${SLURM_CPUS_PER_TASK} --b2-sensitive \
  --GTF=annot/140609_refseq_nomir_nosnor.gtf --no-novel-juncs \    
  /fdb/igenomes/Mus_musculus/UCSC/mm9/Sequence/Bowtie2Index/genome \
  fastq/rnaseq_250k.fastq.gz

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f tophat.swarm -g 10 -t 16 --module tophat
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module tophat Loads the tophat module for each subjob in the swarm 
 | |
 | |
 | |






