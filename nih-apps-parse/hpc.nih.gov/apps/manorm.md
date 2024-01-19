

document.querySelector('title').textContent = 'Manorm on Biowulf';
Manorm on Biowulf


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



 MAnorm is for quantitative comparison of ChIP-Seq data sets describing transcription factor binding sites and epigenetic modifications. The quantitative binding differences inferred by MAnorm showed strong correlation with both the changes in expression of target genes and the binding of cell type-specific regulators.



### References:


* Shao Z, Zhang Y, Yuan GC, Orkin SH, Waxman DJ. (2012) MAnorm: a robust model for quantitative comparison of ChIP-Seq data sets. Genome Biol. Mar 16;13(3):R16.


Documentation
* [Manorm Main Site](https://github.com/shao-lab/MAnorm)


Important Notes
* Module Name: manorm (see [the modules page](/apps/modules.html) for more information)
* singlethreaded app



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem 10g**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load manorm**

[user@cn3144 ~]$ **manorm -h**
usage: manorm [-h] [-v] --p1 PEAKS_FILE1 --p2 PEAKS_FILE2 --r1 READS_FILE1
              --r2 READS_FILE2 [--s1 SHIFT_SIZE1] [--s2 SHIFT_SIZE2]
              [-w WIDTH] [-d DIS_CUTOFF] [-n N_RANDOM] [-m M_CUTOFF]
              [-p P_CUTOFF] [-s] [--name1 NAME1] [--name2 NAME2] -o
              OUTPUT_NAME

MAnorm -- A robust model for quantitative comparison of ChIP-Seq data sets.

optional arguments:
  -h, --help        show this help message and exit
  -v, --version     show program's version number and exit

Input Arguments:
  --p1 PEAKS_FILE1  Peaks file of sample 1. BED and MACS format are currently
                    supported. Please refer to documents for details.
                    (default: None)
  --p2 PEAKS_FILE2  Peaks file of sample 2. (default: None)
  --r1 READS_FILE1  Reads file of sample 1. BED format are currently
                    supported. (default: None)
  --r2 READS_FILE2  Reads file of sample 2. (default: None)
  --s1 SHIFT_SIZE1  Reads shift size of sample 1. This value is used to shift
                    reads towards 3' direction to determine the precise
                    binding site. Set as half of the fragment length.
                    (default: 100)
  --s2 SHIFT_SIZE2  Reads shift size of sample 2. (default: 100)

Normalization Model Arguments:
  -w WIDTH          Width of the window to calculate read densities. Windows
                    with unified length of 2*width centered at peak
                    summit/midpoint are used to quantify the binding signal.
                    This should match the typical length of peaks, a value of
                    1000 is recommended for sharp histone marks like H3K4me3
                    and H3K9/27ac, and 500 for transcription factors or DNase-
                    Seq. (default: 1000)
  -d DIS_CUTOFF     Summit-to-summit distance cutoff for common peaks.
                    Default=width/2. Only overlapped peaks with summit-to-
                    summit distance less than than this value are considered
                    as real common peaks of two samples when fitting M-A
                    normalization model. (default: None)

Advanced Arguments:
  -n N_RANDOM       Number of simulation to test the enrichment of peak
                    overlap between two samples. (default: 10)
  -m M_CUTOFF       M-value cutoff to distinguish biased peaks from unbiased
                    peaks. Peaks with M-value>=M_cutoff and P-value<=P_cutoff
                    are defined as sample1-biased(specific) peaks, while peaks
                    with M-value<=-1*M_cutoff and P-value<=P_cutoff are
                    defined as sample2-biased peaks. (default: 1.0)
  -p P_CUTOFF       P-value cutoff to define biased (sample 1/2-specific)
                    peaks. (default: 0.01)

Output arguments:
  -s                By default, MAnorm will write the comparison results of
                    unique and merged common peaks in a single output file.
                    With this option on, two extra files which contains the
                    results of the original(unmerged) peaks will also be
                    outputted. (default: False)
  --name1 NAME1     Name (experiment condition/cell-type etc.) of sample1. If
                    specified, it will be used to replace the peaks/reads
                    input file name as the sample name in output files.
                    (default: None)
  --name2 NAME2     Name (experiment condition/cell-type etc.) of sample2.
                    (default: None)
  -o OUTPUT_NAME    Output directory. When --name1 and --name2 are not
                    specified, MAnorm will use it as the prefix of comparison
                    output file. (default: None)

[user@cn3144 ~]$ **manorm --p1 sample1\_peaks.bed --p2 sample2\_peaks.bed --r1 sample1\_reads.bed --r2 sample2\_reads.bed -o sample1\_vs\_sample2**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. manorm.sh). For example:



```

#!/bin/bash
set -e
module load manorm
manorm --p1 sample1_peaks.bed --p2 sample2_peaks.bed \
--r1 sample1_reads.bed --r2 sample2_reads.bed \
-o sample1_vs_sample2

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --mem=10g manorm.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. manorm.swarm). For example:



```

manorm --p1 sample1_peaks.bed --p2 sample2_peaks.bed --r1 sample1_reads.bed --r2 sample2_reads.bed -o sample1_vs_sample2
manorm --p1 sample3_peaks.bed --p2 sample4_peaks.bed --r1 sample3_reads.bed --r2 sample4_reads.bed -o sample3_vs_sample4
manorm --p1 sample5_peaks.bed --p2 sample6_peaks.bed --r1 sample5_reads.bed --r2 sample6_reads.bed -o sample5_vs_sample6

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f manorm.swarm -g 10 [-t #] --module manorm
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module manorm Loads the manorm module for each subjob in the swarm 
 | |
 | |
 | |








