

document.querySelector('title').textContent = 'UNICYCLER on Biowulf';
UNICYCLER on Biowulf


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



Unicycler is an assembly pipeline for bacterial genomes. It can assemble Illumina-only read sets where it functions as a SPAdes-optimiser. It can also assembly long-read-only sets (PacBio or Nanopore) where it runs a miniasm+Racon pipeline. 



### References:


* [Unicycler: resolving bacterial genome assemblies from short and long sequencing reads.](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005595) Wick RR, Judd LM, Gorrie CL, Holt KE. PLoS Comput Biol 2017.


Documentation
* [Unicycler Main Site](https://github.com/rrwick/Unicycler)


Important Notes
* Module Name: unicycler (see [the modules page](/apps/modules.html) for more information)



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

[user@cn3144 ~]$ **module load unicycler**
[+] Loading python 3.5  ... 
[+] Loading bowtie  2-2.3.5.1 
[+] Loading racon 1.3.2  ... 
[+] Loading blast 2.9.0+  ... 
[+] Loading pilon  1.23 
[+] Loading unicycler  0.4.8 

[user@cn3144 ~]$ **unicycler --help**
usage: unicycler [-h] [--help_all] [--version] [-1 SHORT1] [-2 SHORT2] [-s UNPAIRED] [-l LONG] -o OUT [--verbosity VERBOSITY] [--min_fasta_length MIN_FASTA_LENGTH]
                 [--keep KEEP] [-t THREADS] [--mode {conservative,normal,bold}] [--linear_seqs LINEAR_SEQS] [--vcf]

       __
       \ \___
        \ ___\
        //
   ____//      _    _         _                     _
 //_  //\\    | |  | |       |_|                   | |
//  \//  \\   | |  | | _ __   _   ___  _   _   ___ | |  ___  _ __
||  (O)  ||   | |  | || '_ \ | | / __|| | | | / __|| | / _ \| '__|
\\    \_ //   | |__| || | | || || (__ | |_| || (__ | ||  __/| |
 \\_____//     \____/ |_| |_||_| \___| \__, | \___||_| \___||_|
                                        __/ |
                                       |___/

Unicycler: an assembly pipeline for bacterial genomes

Help:
  -h, --help                           Show this help message and exit
  --help_all                           Show a help message with all program options
  --version                            Show Unicycler's version number

Input:
  -1 SHORT1, --short1 SHORT1           FASTQ file of first short reads in each pair (required)
  -2 SHORT2, --short2 SHORT2           FASTQ file of second short reads in each pair (required)
  -s UNPAIRED, --unpaired UNPAIRED     FASTQ file of unpaired short reads (optional)
  -l LONG, --long LONG                 FASTQ or FASTA file of long reads (optional)

Output:
  -o OUT, --out OUT                    Output directory (required)
  --verbosity VERBOSITY                Level of stdout and log file information (default: 1)
                                         0 = no stdout, 1 = basic progress indicators, 2 = extra info, 3 = debugging info
  --min_fasta_length MIN_FASTA_LENGTH  Exclude contigs from the FASTA file which are shorter than this length (default: 100)
  --keep KEEP                          Level of file retention (default: 1)
                                         0 = only keep final files: assembly (FASTA, GFA and log), 1 = also save graphs at main checkpoints,
                                         2 = also keep SAM (enables fast rerun in different mode), 3 = keep all temp files and save all graphs (for debugging)
  --vcf                                Produce a VCF by mapping the short reads to the final assembly (experimental, default: do not produce a vcf file)

Other:
  -t THREADS, --threads THREADS        Number of threads used (default: 8)
  --mode {conservative,normal,bold}    Bridging mode (default: normal)
                                         conservative = smaller contigs, lowest misassembly rate
                                         normal = moderate contig size and misassembly rate
                                         bold = longest contigs, higher misassembly rate
  --linear_seqs LINEAR_SEQS            The expected number of linear (i.e. non-circular) sequences in the underlying sequence (default: 0)

[user@cn3144 ~]$ **unicycler -1 short\_reads\_1.fastq.gz -2 short\_reads\_2.fastq.gz -o output\_dir**
Starting Unicycler (2019-11-22 10:42:56)
    Welcome to Unicycler, an assembly pipeline for bacterial genomes. Since you provided only short reads, Unicycler will essentially function as a SPAdes-optimiser. It
will try many k-mer sizes, choose the best based on contig length and graph connectivity, and scaffold the graph using SPAdes repeat resolution.
    For more information, please see https://github.com/rrwick/Unicycler

Command: unicycler -1 short_reads_1.fastq.gz -2 short_reads_2.fastq.gz -o output_dir

Unicycler version: v0.4.8
Using 8 threads

Dependencies:
  Program         Version     Status  
  spades.py       3.13.1      good    
  racon                       not used
  makeblastdb     2.9.0+      good    
  tblastn         2.9.0+      good    
  bowtie2-build   2.3.5.1     good    
  bowtie2         2.3.5.1     good    
  samtools        1.9         good    
  java            1.8.0_181   good    
  pilon           1.23        good    
  bcftools                    not used


SPAdes read error correction (2019-11-22 10:42:59)
    Unicycler uses the SPAdes read error correction module to reduce the number of errors in the short read before SPAdes assembly. This can make the assembly faster and
simplify the assembly graph structure.

Command: spades.py -1 TEST_DATA/short_reads_1.fastq.gz -2 TEST_DATA/short_reads_2.fastq.gz -o TEST_DATA/output_dir/spades_assembly/read_correction --threads 8 --only-error-correction

Corrected reads:
  TEST_DATA/output_dir/spades_assembly/corrected_1.fastq.gz
  TEST_DATA/output_dir/spades_assembly/corrected_2.fastq.gz


Choosing k-mer range for assembly (2019-11-22 10:43:38)
    Unicycler chooses a k-mer range for SPAdes based on the length of the input reads. It uses a wide range of many k-mer sizes to maximise the chance of finding an ideal
assembly.

SPAdes maximum k-mer: 127
Median read length: 125
K-mer range: 25, 43, 59, 73, 83, 93, 101, 107, 113, 119


SPAdes assemblies (2019-11-22 10:43:39)
    Unicycler now uses SPAdes to assemble the short reads. It scores the assembly graph for each k-mer using the number of contigs (fewer is better) and the number of dead
ends (fewer is better). The score function is 1/(c*(d+2)), where c is the contig count and d is the dead end count.

K-mer   Contigs   Dead ends   Score   
   25                           failed
   43                           failed
   59                           failed
   73                           failed
   83                           failed
   93                           failed
  101                           failed
  107                           failed
  113                           failed
  119       111          28   3.00e-04 ← best

Read depth filter: removed 0 contigs totalling 0 bp
Deleting TEST_DATA/output_dir/spades_assembly/


Determining graph multiplicity (2019-11-22 10:44:38)
    Multiplicity is the number of times a sequence occurs in the underlying sequence. Single-copy contigs (those with a multiplicity of one, occurring only once in the
underlying sequence) are particularly useful.

Saving TEST_DATA/output_dir/001_best_spades_graph.gfa


Cleaning graph (2019-11-22 10:44:38)
    Unicycler now performs various cleaning procedures on the graph to remove overlaps and simplify the graph structure. The end result is a graph ready for bridging.

Graph overlaps removed

Removed zero-length segments:
    83, 84, 86, 87, 91, 94, 100, 103

Removed zero-length segments:
    85

Merged small segments:
    108, 109, 110, 111

Saving TEST_DATA/output_dir/002_overlaps_removed.gfa

    Unicycler now selects a set of anchor contigs from the single-copy contigs. These are the contigs which will be connected via bridges to form the final assembly.

36 anchor segments (174,096 bp) out of 98 total segments (191,318 bp)


Creating SPAdes contig bridges (2019-11-22 10:44:38)
    SPAdes uses paired-end information to perform repeat resolution (RR) and produce contigs from the assembly graph. SPAdes saves the graph paths corresponding to these
contigs in the contigs.paths file. When one of these paths contains two or more anchor contigs, Unicycler can create a bridge from the path.

                                 Bridge
Start        Path       End     quality
   -4                   12         55.1
    1                   -22        62.8
    9        -34        1           9.2
   11      77 → -62     22         59.9
   12         47        28         54.7
   18   -77 → 76 → 80   33         42.8
   31        -72        -32        62.2
   33      -80 → 85     30         51.2
   37         75        11         60.3


Creating loop unrolling bridges (2019-11-22 10:44:38)
    When a SPAdes contig path connects an anchor contig with the middle contig of a simple loop, Unicycler concludes that the sequences are contiguous (i.e. the loop is
not a separate piece of DNA). It then uses the read depth of the middle and repeat contigs to guess the number of times to traverse the loop and makes a bridge.

No loop unrolling bridges made


Applying bridges (2019-11-22 10:44:38)
    Unicycler now applies to the graph in decreasing order of quality. This ensures that when multiple, contradictory bridges exist, the most supported option is used.

Bridge type   Start → end   Path          Quality
SPAdes            1 → -22                  62.822
SPAdes           31 → -32   -72            62.196
SPAdes           37 → 11    75             60.300
SPAdes           11 → 22    77, -62        59.937
SPAdes           -4 → 12                   55.117
SPAdes           12 → 28    47             54.728
SPAdes           33 → 30    -80, 85        51.155
SPAdes           18 → 33    -77, 76, 80    42.846

Saving TEST_DATA/output_dir/003_bridges_applied.gfa


Bridged assembly graph (2019-11-22 10:44:38)
    The assembly is now mostly finished and no more structural changes will be made. Ideally the assembly graph should now have one contig per replicon and no erroneous
contigs (i.e a complete assembly). If there are more contigs, then the assembly is not complete.

Saving TEST_DATA/output_dir/004_final_clean.gfa

Component   Segments   Links   Length    N50     Longest segment   Status    
    total         66      77   189,317   7,866            43,777             
        1         63      76   178,782   8,115            43,777   incomplete
        2          1       1     5,153   5,153             5,153     complete
        3          1       0     3,283   3,283             3,283   incomplete
        4          1       0     2,099   2,099             2,099   incomplete


Polishing assembly with Pilon (2019-11-22 10:44:38)
    Unicycler now conducts multiple rounds of Pilon in an attempt to repair any remaining small-scale errors with the assembly.

Aligning reads to find appropriate insert size range...
Insert size 1st percentile:  264
Insert size 99th percentile: 549

Pilon polish round 1
Unable to polish assembly using Pilon: Pilon did not produce FASTA file


Rotating completed replicons (2019-11-22 10:45:09)
    Any completed circular contigs (i.e. single contigs which have one link connecting end to start) can have their start position changed without altering the sequence.
For consistency, Unicycler now searches for a starting gene (dnaA or repA) in each such contig, and if one is found, the contig is rotated to start with that gene on the
forward strand.

Segment   Length   Depth   Starting gene   Position   Strand   Identity   Coverage
     11    5,153   4.01x   none found                                             


Assembly complete (2019-11-22 10:45:23)
Saving TEST_DATA/output_dir/assembly.gfa
Saving TEST_DATA/output_dir/assembly.fasta

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. unicycler.sh). For example:



```

#!/bin/bash
# sbatch --gres=lscratch:10 --cpus-per-task=16 unicycler.sh 
set -e
module load unicycler
cp /data/$SLURM_JOB_USER/UNICYCLER_TEST/TEST_DATA/short_reads_1.fastq.gz /lscratch/$SLURM_JOB_ID/.
cp /data/$SLURM_JOB_USER/UNICYCLER_TEST/TEST_DATA/short_reads_2.fastq.gz /lscratch/$SLURM_JOB_ID/.
unicycler -1 /lscratch/$SLURM_JOB_ID/short_reads_1.fastq.gz \
          -2 /lscratch/$SLURM_JOB_ID/short_reads_2.fastq.gz \
          -o /lscratch/$SLURM_JOB_ID/OUTPUT
cp -fr /lscratch/$SLURM_JOB_ID/OUTPUT /data/$SLURM_JOB_USER/UNICYCLER_TEST/TEST_DATA/.

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--gres=lscratch:#] [--cpus-per-task=#] [--mem=#] unicycler.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. unicycler.swarm). For example:



```

#
# swarm -f unicycler.swarm --gres=lscratch:10 --threads-per-process=16 --module unicycler
#
cp /data/$SLURM_JOB_USER/UNICYCLER_TEST/TEST_DATA/short_reads_1.fastq.gz /lscratch/$SLURM_JOB_ID/; \
cp /data/$SLURM_JOB_USER/UNICYCLER_TEST/TEST_DATA/short_reads_2.fastq.gz /lscratch/$SLURM_JOB_ID/; \
unicycler -1 /lscratch/$SLURM_JOB_ID/short_reads_1.fastq.gz \
          -2 /lscratch/$SLURM_JOB_ID/short_reads_2.fastq.gz \
          -o /lscratch/$SLURM_JOB_ID/OUTPUT; \
mv /lscratch/$SLURM_JOB_ID/OUTPUT /data/$SLURM_JOB_USER/UNICYCLER_TEST/TEST_DATA/.
cp /data/$SLURM_JOB_USER/UNICYCLER_TEST/TEST_DATA/long_reads_high_depth.fastq.gz /lscratch/$SLURM_JOB_ID/ ; \
unicycler -l /lscratch/$SLURM_JOB_ID/long_reads_high_depth.fastq.gz \
          -o /lscratch/$SLURM_JOB_ID/OUTPUT2; \
mv /lscratch/$SLURM_JOB_ID/OUTPUT2 /data/$SLURM_JOB_USER/UNICYCLER_TEST/TEST_DATA/.

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f unicycler.swarm [--gres=lscratch:#] [-g #] [-t #] --module unicycler
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module unicycler Loads the unicycler module for each subjob in the swarm 
 | |
 | |
 | |








