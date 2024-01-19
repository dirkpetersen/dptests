

document.querySelector('title').textContent = 'seqkit on Biowulf';
seqkit on Biowulf


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



Seqkit is a rapid tool for manipulating fasta and fastq files. It includes a number of different tools:
format conversion, searching, bam processing and monitoring, filtering and ordering. SeqKit demonstrates competitive 
performance in execution time and memory usage compared to similar tools.



Documentation
* seqkit on [Github](https://github.com/shenwei356/seqkit)
* [Documentation](https://bioinf.shenwei.me/seqkit/)


Important Notes
* Module Name: seqkit (see [the modules page](/apps/modules.html) for more information)
* Example files in `$SEQKIT_TEST_DATA`
* Some fasta/fastq manipulations can be parallelized with `-j/--threads`. Please match
 the number of allocated CPUs to the number of threads.
* Seqkit's memory occupation is small, because most subcommands do not read whole FASTA/Q records in to memory.



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

[user@cn3144]$ **module load seqkit**
[user@cn3144]$ **seqkit help** 
seqKit -- a cross-platform and ultrafast toolkit for FASTA/Q file manipulation

Version: 0.12.1

Author: Wei Shen 

Documents : http://bioinf.shenwei.me/seqkit
Source code: https://github.com/shenwei356/seqkit
Please cite: https://doi.org/10.1371/journal.pone.0163962

Usage:
 seqkit [command]

Available Commands:
 amplicon retrieve amplicon (or specific region around it) via primer(s)
 bam monitoring and online histograms of BAM record features
 common find common sequences of multiple files by id/name/sequence
 concat concatenate sequences with same ID from multiple files
 convert convert FASTQ quality encoding between Sanger, Solexa and Illumina
 duplicate duplicate sequences N times
 faidx create FASTA index file and extract subsequence
 fish look for short sequences in larger sequences using local alignment
 fq2fa convert FASTQ to FASTA
 fx2tab convert FASTA/Q to tabular format (with length/GC content/GC skew)
 genautocomplete generate shell autocompletion script
 grep search sequences by ID/name/sequence/sequence motifs, mismatch allowed
 head print first N FASTA/Q records
 help Help about any command
 locate locate subsequences/motifs, mismatch allowed
 mutate edit sequence (point mutation, insertion, deletion)
 range print FASTA/Q records in a range (start:end)
 rename rename duplicated IDs
 replace replace name/sequence by regular expression
 restart reset start position for circular genome
 rmdup remove duplicated sequences by id/name/sequence
 sample sample sequences by number or proportion
 sana sanitize broken single line fastq files
 seq transform sequences (revserse, complement, extract ID...)
 shuffle shuffle sequences
 sliding sliding sequences, circular genome supported
 sort sort sequences by id/name/sequence/length
 split split sequences into files by id/seq region/size/parts (mainly for FASTA)
 split2 split sequences into files by size/parts (FASTA, PE/SE FASTQ)
 stats simple statistics of FASTA/Q files
 subseq get subsequences by region/gtf/bed, including flanking sequences
 tab2fx convert tabular format to FASTA/Q format
 translate translate DNA/RNA to protein sequence (supporting ambiguous bases)
 version print version information and check for update
 watch monitoring and online histograms of sequence features

Flags:
 --alphabet-guess-seq-length int length of sequence prefix of the first FASTA record based on which seqkit guesses the sequence type (0 for whole seq) (default 10000)
 -h, --help help for seqkit
 --id-ncbi FASTA head is NCBI-style, e.g. >gi|110645304|ref|NC\_002516.2| Pseud...
 --id-regexp string regular expression for parsing ID (default "^(\\S+)\\s?")
 --infile-list string file of input files list (one file per line), if given, they are appended to files from cli arguments
 -w, --line-width int line width when outputing FASTA format (0 for no wrap) (default 60)
 -o, --out-file string out file ("-" for stdout, suffix .gz for gzipped out) (default "-")
 --quiet be quiet and do not show extra information
 -t, --seq-type string sequence type (dna|rna|protein|unlimit|auto) (for auto, it automatically detect by the first sequence) (default "auto")
 -j, --threads int number of CPUs. (default value: 1 for single-CPU PC, 2 for others) (default 2)

Use "seqkit [command] --help" for more information about a command.
[user@cn3144]$ **cp $SEQKIT\_TEST\_DATA/\* .**
[user@cn3144]$ **ls -lh**
total 13M
-rw-rw-r-- 1 user group 1.5M Mar 11 2018 hairpin.fa.gz
-rw-rw-r-- 1 user group 12M Sep 10 14:51 read1\_250k.fastq.gz

[user@cn3144]$ **seqkit stat hairpin.fa.gz**
file format type num\_seqs sum\_len min\_len avg\_len max\_len
hairpin.fa.gz FASTA RNA 38,589 3,729,811 39 96.7 2,354
[user@cn3144]$ **zcat hairpin.fa.gz|seqkit grep -s -i -p AAUCUUCCUUUGUCU -m 1**
>mdo-mir-12331 MI0040833 Monodelphis domestica miR-12331 stem-loop
AAUCUUCCUUUGUCUACUUUAAGUUCCUGAUUGACUAACUUAAAAUAGACAAAGUGAGAG
UU
[user@cn3144]$ **gunzip hairpin.fa.gz**
[user@cn3144]$ **seqkit sort --by-length --reverse --two-pass hairpin.fa > hairpin\_sorted.fa**
[user@cn3144]$ **seqkit fx2tab --length hairpin\_sorted.fa --name --only-id | head**
atr-MIR8591 2354
eun-MIR10219 1545
atr-MIR8598 1411
aly-MIR858 938
mdm-MIR858 921
mtr-MIR7700 910
cre-MIR914 894
eun-MIR10211a 820
eun-MIR10211b 820
atr-MIR8612 805

[user@cn3144]$ **seqkit rmdup -s -i hairpin.fa -o clean.fa.gz -D dups.txt**
[user@cn3144]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. seqkit.sh), which uses the input file 'seqkit.in'. For example:



```

#!/bin/bash
module load seqkit || exit 1
seqkit split $SEQKIR_TEST_DATA/hairpin.fa.gz -i --id-regexp "^([\w]+)\-" --two-pass

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=2 --mem=4g seqkit.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. seqkit.swarm). For example:



```

seqkit -B -g -G -l -n file1.fa
seqkit -B -g -G -l -n file2.fa
seqkit -B -g -G -l -n file3.fa

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f seqkit.swarm -g 4 -t 1 --module seqkit
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g #  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t #  Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module seqkit  Loads the seqkit module for each subjob in the swarm 
 | |
 | |
 | |








