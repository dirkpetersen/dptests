

document.querySelector('title').textContent = 'DETONATE: DE novo TranscriptOme rNa-seq Assembly with or without the Truth Evaluation ';
**DETONATE: DE novo TranscriptOme rNa-seq Assembly with or without the Truth Evaluation** 


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



DETONATE consists of two component packages for evaluation of de novo transcriptome assemblies, RSEM-EVAL and REF-EVAL. RSEM-EVAL is a reference-free evaluation method based on a novel probabilistic model that depends only on an assembly and the RNA-Seq reads used for its construction. REF-EVAL is a toolkit of reference-based measures. 



### References:


* Bo Li, Nathanael Fillmore, Yongsheng Bai, Mike Collins, James A. Thomson, Ron Stewart, and Colin N. Dewey.   

 *Evaluation of de novo transcriptome assemblies from RNA-Seq data*  

[Genome Biology](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0553-5)  2014, **15**:553


Documentation
* [DETONATE Home page](http://deweylab.biostat.wisc.edu/detonate/)


Important Notes
* Module Name: DETONATE (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set
	+ **DETONATE\_HOME**  installation directory
	+ **DETONATE\_BIN**       executable directory
	+ **DETONATE\_SRC**       source code directory
	+ **DETONATE\_DATA**  sample data directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=4g**
[user@@cn3200 ~]$**module load DETONATE** 
[+] Loading gcc  7.2.0  ... 
[+] Loading bowtie  1.2.2 
[+] Loading GSL 2.4 for GCC 7.2.0 ... 
[+] Loading openmpi 3.0.0  for GCC 7.2.0 
[+] Loading R 3.5.0_build2
[+] Loading DETONATE  1.11

```

Available executables are: 

```

**ref-eval
ref-eval-estimate-true-assembly
rsem-build-read-index
rsem-eval-calculate-score
rsem-eval-estimate-transcript-length-distribution
rsem-eval-run-em
rsem-extract-reference-transcripts
rsem-parse-alignments
rsem-plot-model
rsem-preref
rsem-sam-validator
rsem-scan-for-paired-end-reads
rsem-simulate-reads
rsem-synthesis-reference-transcripts**    

```

In order to display a usage message for an executable, type its name.   

For eample:

```

[user@biowulf]$ **ref-eval**
  REF-EVAL: A toolkit of reference-based scores for de novo transcriptome
                       sequence assembly evaluation

Overview

   REF-EVAL computes a number of reference-based scores. These scores
   measure the quality of a transcriptome assembly relative to a
   collection of reference sequences. For information about how to run
   REF-EVAL, see "Usage" and the sections following it below. For
   information about the score definitions, see "Score definitions" and
   the sections following it below.

Usage

   As an optional first step, estimate the "true" assembly, using
   [1]REF-EVAL-ESTIMATE-TRUE-ASSEMBLY. Alternatively, you can use the
   full-length reference transcript sequences directly as a reference.

   From now on, we will call the estimated "true" assembly or the
   collection of full-length reference sequences (whichever you choose
   to use) the reference. Let's assume that the assembly of interest is
   in A.fa, and the reference is in B.fa.
...
 $ ./ref-eval --scores=nucl,pair,contig,kmer,kc \
              --weighted=both \
              --A-seqs A.fa \
              --B-seqs B.fa \
              --A-expr A_expr.isoforms.results \
              --B-expr B_expr.isoforms.results \
              --A-to-B A_to_B.psl \
              --B-to-A B_to_A.psl \
              --num-reads 5000000 \
              --readlen 76 \
              --kmerlen 76 \
              | tee scores.txt
...
[user@@cn3200 ~]$ **rsem-eval-calculate-score** 
SYNOPSIS
     rsem-eval-calculate-score [options] upstream_read_file(s) assembly_fasta_file sample_name L
     rsem-eval-calculate-score [options] --paired-end upstream_read_file(s) downstream_read_file(s) assembly_fasta_file sample_name L
     rsem-eval-calculate-score [options] --sam/--bam [--paired-end] input assembly_fasta_file sample_name L

ARGUMENTS
...

BASIC OPTIONS
...

```

Here is a sample command for running executable rsem-eval-calculate-score:

```

[user@@cn3200 ~]$ **rsem-eval-calculate-score $DETONATE\_DATA/toy\_SE.fq $DETONATE\_DATA/toy\_assembly\_1.fa ./rsem\_eval\_1 76 --transcript-length-parameters $DETONATE\_HOME/src/true\_transcript\_length\_distribution/mouse.txt -p 16**
me

rsem-synthesis-reference-transcripts ./rsem_eval_1.temp/rsem_eval_1 0 0 0 /usr/local/apps/DETONATE/1.11/sample_data/toy_assembly_1.fa
Transcript Information File is generated!
Group File is generated!
Extracted Sequences File is generated!

rsem-preref ./rsem_eval_1.temp/rsem_eval_1.transcripts.fa 1 ./rsem_eval_1.temp/rsem_eval_1
Refs.makeRefs finished!
Refs.saveRefs finished!
./rsem_eval_1.temp/rsem_eval_1.idx.fa is generated!
./rsem_eval_1.temp/rsem_eval_1.n2g.idx.fa is generated!

bowtie-build -f ./rsem_eval_1.temp/rsem_eval_1.n2g.idx.fa ./rsem_eval_1.temp/rsem_eval_1
Settings:
  Output files: "./rsem_eval_1.temp/rsem_eval_1.*.ebwt"
  Line rate: 6 (line is 64 bytes)
  Lines per side: 1 (side is 64 bytes)
  Offset rate: 5 (one in 32)
  FTable chars: 10
  Strings: unpacked
  Max bucket size: default
  Max bucket size, sqrt multiplier: default
  Max bucket size, len divisor: 4
  Difference-cover sample period: 1024
  Endianness: little
  Actual local endianness: little
  Sanity checking: disabled
  Assertions: disabled
  Random seed: 0
  Sizeofs: void*:8, int:4, long:8, size_t:8
Input files DNA, FASTA:
  ./rsem_eval_1.temp/rsem_eval_1.n2g.idx.fa
Reading reference sizes
  Time reading reference sizes: 00:00:00
Calculating joined length
Writing header
Reserving space for joined string
Joining reference sequences
  Time to join reference sequences: 00:00:00
bmax according to bmaxDivN setting: 156
Using parameters --bmax 117 --dcv 1024
  Doing ahead-of-time memory usage test
  Passed!  Constructing with these parameters: --bmax 117 --dcv 1024
Constructing suffix-array element generator
Building DifferenceCoverSample
  Building sPrime
  Building sPrimeOrder
  V-Sorting samples
  V-Sorting samples time: 00:00:00
  Allocating rank array
  Ranking v-sort output
  Ranking v-sort output time: 00:00:00
  Invoking Larsson-Sadakane on ranks
  Invoking Larsson-Sadakane on ranks time: 00:00:00
  Sanity-checking and returning
Building samples
Reserving space for 12 sample suffixes
Generating random suffixes
QSorting 12 sample offsets, eliminating duplicates
QSorting sample offsets, eliminating duplicates time: 00:00:00
Multikey QSorting 11 samples
  (Using difference cover)
  Multikey QSorting samples time: 00:00:00
Calculating bucket sizes
Splitting and merging
  Splitting and merging time: 00:00:00
Avg bucket size: 626 (target: 116)
Converting suffix-array elements to index image
Allocating ftab, absorbFtab
Entering Ebwt loop
Getting block 1 of 1
  No samples; assembling all-inclusive block
  Sorting block of length 626 for bucket 1
  (Using difference cover)
  Sorting block time: 00:00:00
Returning block of 627 for bucket 1
Exited Ebwt loop
fchr[A]: 0
fchr[C]: 190
fchr[G]: 340
fchr[T]: 501
fchr[$]: 626
Exiting Ebwt::buildToDisk()
Returning from initFromVector
Wrote 4194735 bytes to primary EBWT file: ./rsem_eval_1.temp/rsem_eval_1.1.ebwt
Wrote 84 bytes to secondary EBWT file: ./rsem_eval_1.temp/rsem_eval_1.2.ebwt
Re-opening _in1 and _in2 as input streams
Returning from Ebwt constructor
Headers:
    len: 626
    bwtLen: 627
    sz: 157
    bwtSz: 157
    lineRate: 6
    linesPerSide: 1
    offRate: 5
    offMask: 0xffffffe0
    isaRate: -1
    isaMask: 0xffffffff
    ftabChars: 10
    eftabLen: 20
    eftabSz: 80
    ftabLen: 1048577
    ftabSz: 4194308
    offsLen: 20
    offsSz: 80
    isaLen: 0
    isaSz: 0
    lineSz: 64
    sideSz: 64
    sideBwtSz: 56
    sideBwtLen: 224
    numSidePairs: 2
    numSides: 4
    numLines: 4
    ebwtTotLen: 256
    ebwtTotSz: 256
    reverse: 0
Total time for call to driver() for forward index: 00:00:00
Reading reference sizes
  Time reading reference sizes: 00:00:00
Calculating joined length
Writing header
Reserving space for joined string
Joining reference sequences
  Time to join reference sequences: 00:00:00
bmax according to bmaxDivN setting: 156
Using parameters --bmax 117 --dcv 1024
  Doing ahead-of-time memory usage test
  Passed!  Constructing with these parameters: --bmax 117 --dcv 1024
Constructing suffix-array element generator
Building DifferenceCoverSample
  Building sPrime
  Building sPrimeOrder
  V-Sorting samples
  V-Sorting samples time: 00:00:00
  Allocating rank array
  Ranking v-sort output
  Ranking v-sort output time: 00:00:00
  Invoking Larsson-Sadakane on ranks
  Invoking Larsson-Sadakane on ranks time: 00:00:00
  Sanity-checking and returning
Building samples
Reserving space for 12 sample suffixes
Generating random suffixes
QSorting 12 sample offsets, eliminating duplicates
QSorting sample offsets, eliminating duplicates time: 00:00:00
Multikey QSorting 11 samples
  (Using difference cover)
  Multikey QSorting samples time: 00:00:00
Calculating bucket sizes
Splitting and merging
  Splitting and merging time: 00:00:00
Avg bucket size: 626 (target: 116)
Converting suffix-array elements to index image
Allocating ftab, absorbFtab
Entering Ebwt loop
Getting block 1 of 1
  No samples; assembling all-inclusive block
  Sorting block of length 626 for bucket 1
  (Using difference cover)
  Sorting block time: 00:00:00
Returning block of 627 for bucket 1
Exited Ebwt loop
fchr[A]: 0
fchr[C]: 190
fchr[G]: 340
fchr[T]: 501
fchr[$]: 626
Exiting Ebwt::buildToDisk()
Returning from initFromVector
Wrote 4194735 bytes to primary EBWT file: ./rsem_eval_1.temp/rsem_eval_1.rev.1.ebwt
Wrote 84 bytes to secondary EBWT file: ./rsem_eval_1.temp/rsem_eval_1.rev.2.ebwt
Re-opening _in1 and _in2 as input streams
Returning from Ebwt constructor
Headers:
    len: 626
    bwtLen: 627
    sz: 157
    bwtSz: 157
    lineRate: 6
    linesPerSide: 1
    offRate: 5
    offMask: 0xffffffe0
    isaRate: -1
    isaMask: 0xffffffff
    ftabChars: 10
    eftabLen: 20
    eftabSz: 80
    ftabLen: 1048577
    ftabSz: 4194308
    offsLen: 20
    offsSz: 80
    isaLen: 0
    isaSz: 0
    lineSz: 64
    sideSz: 64
    sideBwtSz: 56
    sideBwtLen: 224
    numSidePairs: 2
    numSides: 4
    numLines: 4
    ebwtTotLen: 256
    ebwtTotSz: 256
    reverse: 0
Total time for backward call to driver() for mirror index: 00:00:01

bowtie -q --phred33-quals -n 2 -e 99999999 -l 25 -p 16 -a -m 200 -S ./rsem_eval_1.temp/rsem_eval_1 /usr/local/apps/DETONATE/1.11/sample_data/toy_SE.fq | samtools view -S -b -o ./rsem_eval_1.temp/rsem_eval_1.bam -
[samopen] SAM header is present: 1 sequences.
# reads processed: 3818
# reads with at least one reported alignment: 3812 (99.84%)
# reads that failed to align: 6 (0.16%)
Reported 3812 alignments

rsem-parse-alignments ./rsem_eval_1.temp/rsem_eval_1 ./rsem_eval_1.temp/rsem_eval_1 ./rsem_eval_1.stat/rsem_eval_1 b ./rsem_eval_1.temp/rsem_eval_1.bam -t 1 -tag XM
Done!

rsem-build-read-index 32 1 0 ./rsem_eval_1.temp/rsem_eval_1_alignable.fq
Build Index ./rsem_eval_1.temp/rsem_eval_1_alignable.fq is Done!

rsem-eval-run-em ./rsem_eval_1.temp/rsem_eval_1 1 ./rsem_eval_1 ./rsem_eval_1.temp/rsem_eval_1 ./rsem_eval_1.stat/rsem_eval_1 1.00636237352512811 0.00047971891527351895 76 0 -p 16
Refs.loadRefs finished!
Thread 0 : N = 238, NHit = 238
Thread 1 : N = 238, NHit = 238
Thread 2 : N = 238, NHit = 238
Thread 3 : N = 238, NHit = 238
Thread 4 : N = 238, NHit = 238
Thread 5 : N = 238, NHit = 238
Thread 6 : N = 238, NHit = 238
Thread 7 : N = 238, NHit = 238
Thread 8 : N = 238, NHit = 238
Thread 9 : N = 238, NHit = 238
Thread 10 : N = 238, NHit = 238
Thread 11 : N = 238, NHit = 238
Thread 12 : N = 238, NHit = 238
Thread 13 : N = 238, NHit = 238
Thread 14 : N = 238, NHit = 238
DAT 0 reads left
Thread 15 : N = 242, NHit = 242
EM_init finished!
estimateFromReads, N0 finished.
estimateFromReads, N1 finished.
ROUND = 1, SUM = 3818, bChange = 2.52329e-05, totNum = 0
ROUND = 2, SUM = 3818, bChange = 2.51884e-05, totNum = 0
ROUND = 3, SUM = 3818, bChange = 1.5407e-11, totNum = 0
ROUND = 4, SUM = 3818, bChange = 0, totNum = 0
ROUND = 5, SUM = 3818, bChange = 0, totNum = 0
ROUND = 6, SUM = 3818, bChange = 0, totNum = 0
ROUND = 7, SUM = 3818, bChange = 0, totNum = 0
ROUND = 8, SUM = 3818, bChange = 0, totNum = 0
ROUND = 9, SUM = 3818, bChange = 0, totNum = 0
ROUND = 10, SUM = 3818, bChange = 0, totNum = 0
ROUND = 11, SUM = 3818, bChange = 0, totNum = 0
ROUND = 12, SUM = 3818, bChange = 0, totNum = 0
ROUND = 13, SUM = 3818, bChange = 0, totNum = 0
ROUND = 14, SUM = 3818, bChange = 0, totNum = 0
ROUND = 15, SUM = 3818, bChange = 0, totNum = 0
ROUND = 16, SUM = 3818, bChange = 0, totNum = 0
ROUND = 17, SUM = 3818, bChange = 0, totNum = 0
ROUND = 18, SUM = 3818, bChange = 0, totNum = 0
ROUND = 19, SUM = 3818, bChange = 0, totNum = 0
ROUND = 20, SUM = 3818, bChange = 0, totNum = 0
Expression Results are written!
ces is created.
Calculating assembly priors.
Assembly priors are calculated!
Calculating maximum data likelihood in log space.
Maximum data likelihood is calculated!
Calculating data likelihood in log space.
Data likelihood is calculated!
Calculating correction score.
Correction score is calculated!
score is written.
Time Used for EM.cpp : 0 h 00 m 00 s

rm -rf ./rsem_eval_1.temp



```

End the interactive session:

```

[user@cn3200 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```

Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. detonate.sh). For example:



```

#!/bin/bash
#SBATCH --mem=4g
module load DETONATE             
cd /data/$USER
rsem-eval-calculate-score $DETONATE_DATA/toy_SE.fq $DETONATE_DATA/toy_assembly_1.fa ./rsem_eval_1 76 --transcript-length-parameters $DETONATE_DATA/true_transcript_length_distribution/mouse.txt -p 16

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch detonate.sh 
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. detonate.swarm). For example:



```

#!/bin/bash
module load DETONATE
cd /data/$USER
rsem-eval-calculate-score $DETONATE_DATA/toy_SE.fq $DETONATE_DATA/toy_assembly_1.fa ./rsem_eval_1 76 --transcript-length-parameters $DETONATE_DATA/true_transcript_length_distribution/mouse.txt -p 16

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f detonate.swarm -g 4
```





