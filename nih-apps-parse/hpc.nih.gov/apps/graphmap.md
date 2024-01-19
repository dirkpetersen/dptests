

document.querySelector('title').textContent = 'GraphMap - A highly sensitive and accurate mapper for long, error-prone reads';
**GraphMap - A highly sensitive and accurate mapper for long, error-prone reads**


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



GraphMap is a highly sensitive and accurate mapper for long, error-prone reads. 
It offers a number of valuable features, such as mapping position agnostic to alignment parameters,
high sensitivity and precision, handling circular genomes, meaningful mapping quality,
various alignment strategies, and more.



### References:


* Ivan Sović, Mile Šikić, Andreas Wilm, Shannon Nicole Fenlon, Swaine Chen & Niranjan Nagarajan   

*Fast and sensitive mapping of nanopore sequencing reads with GraphMap,*   

[Nature Communications volume 7, Article number: 11307 (2016)](https://www.nature.com/articles/ncomms11307) .


Documentation
* [GraphMap GitHub page](https://github.com/marbl/graphmap)


Important Notes
* Module Name: graphmap (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Implemented as a Singularity container
* Unusual environment variables set
	+ **GM\_HOME**  installation directory
	+ **GM\_HOME**  executable directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=4g**
[user@cn3316 ~]$ **module load graphmap** 
[+] Loading graphmap 0.5.2  ...
[user@cn3316 ~]$ **graphmap -h**
...

Usage:
  graphmap tool

Options
    tool       STR   Specifies the tool to run:
                       align - the entire GraphMap pipeline.
                       owler - Overlapping With Long Erroneous Reads.

```

The grapgmap application can be used as follows:   
   

\*\*Align\*\* all reads from a given FASTA/FASTQ file using anchored alignment approach:

```

[user@cn3316 ~]$  **graphmap align -r escherichia\_coli.fa -d reads.fastq -o alignments.sam** 

```

\*\*Overlap\*\* all reads from a given FASTA/FASTQ file and report overlaps in MHAP format (fast): 

```

[user@cn3316 ~]$ **graphmap owler -r reads.fa -d reads.fa -o overlaps.mhap** 

```

\*\*Align\*\* all reads to a transcriptome sequence: 

```

[user@cn3316 ~]$ **graphmap align -r scerevisiae.fa --gtf scerevisiae.gtf -d reads.fastq -o alignments.sam**  

```

Align all reads and report alignments using the extended CIGAR format. 

```

[user@cn3316 ~]$ **graphmap align -r escherichia\_coli.fa -d reads.fastq -o alignments.sam --extcigar**  

```

Align all reads from a given FASTA/FASTQ file with default number of threads using semiglobal bit-vector alignment: 

```

[user@cn3316 ~]$ **graphmap align -a sg -r escherichia\_coli.fa -d reads.fastq -o alignments.sam** 

```

Overlap all reads from a given FASTA/FASTQ in a full GraphMap mode with generating alignments (slow): 

```

[user@cn3316 ~]$ **graphmap align -x overlap -r reads.fa -d reads.fa -o overlaps.sam**  

```

Align reads using the Gotoh for semiglobal alignment: 

```

[user@cn3316 ~]$ **graphmap align -a sggotoh -r escherichia\_coli.fa -d reads.fastq -o alignments.sam** 

```

Align reads using Gotoh alignment with anchored approach: 

```

[user@cn3316 ~]$ **graphmap align -a anchorgotoh -r escherichia\_coli.fa -d reads.fastq -o alignments.sam** 

```

Process reads from a circular genome: 

```

[user@cn3316 ~]$ **graphmap align -C -r escherichia\_coli.fa -d reads.fastq -o alignments.sam** 

```

Threshold the E-value of alignments to 1e-100. Alignments with E-value > 1e-100 will be called unmapped: 

```

[user@cn3316 ~]$ **graphmap align --evalue 1e-100 -r escherichia\_coli.fa -d reads.fastq -o alignments.sam**  

```

Output all secondary alignments instead of only one best: 

```

[user@cn3316 ~]$ **graphmap align --secondary -r escherichia\_coli.fa -d reads.fastq -o alignments.sam** 

```

Control the similarity for secondary alignments. All alignments to within F\*num\_covered\_bases from the best will be output. 

```

[user@cn3316 ~]$ **graphmap align --secondary -F 0.05 -r escherichia\_coli.fa -d reads.fastq -o alignments.sam** 

```

Limit the number of threads to 8, and load reads in batches of 50MB: 

```

[user@cn3316 ~]$ **graphmap align -t 8 -B 50 -r escherichia\_coli.fa -d reads.fastq -o alignments.sam**  

```

Align reads using more sensitive parameters for Illumina data: 

```

[user@cn3316 ~]$ **graphmap align -x illumina -r escherichia\_coli.fa -d reads.fastq -o alignments.sam** 

```

Load all reads in one batch and align only the first 1000 reads: 

```

[user@cn3316 ~]$ **graphmap align -B 0 -n 1000 -r escherichia\_coli.fa -d reads.fastq -o alignments.sam**  

```

Rebuild the index if it already exists: 

```

[user@cn3316 ~]$ **graphmap align --rebuild-index -r escherichia\_coli.fa -d reads.fastq -o alignments.sam**  

```

Generate only the index. 

```

[user@cn3316 ~]$ **graphmap align -I -r escherichia\_coli.fa** 

```

Run a debug version of GraphMap (build with "make debug") and verbose the SAM output to see various info about alignment: 

```

[user@cn3316 ~]$ **graphmap-debug align -b 3 -r escherichia\_coli.fa -d reads.fastq -o alignments.sam** 

```


End the interactive session:

```

[user@cn3316 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```





