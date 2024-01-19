

document.querySelector('title').textContent = 'bedtools on Biowulf';
bedtools on Biowulf


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


The BEDTools utilities allow one to address common genomics tasks such
finding feature overlaps and computing coverage. In addition, one can develop
sophisticated pipelines that answer complicated research questions by
chaining several BEDTools together.


Bedtools contains a single binary (`bedtools`) that takes a 
number of commands as as it's first argument.



```

bedtools <subcommand> [options]

```

###  Subcommands available in version 2.29.0





| **Genome** |
| [intersect](http://bedtools.readthedocs.org/en/latest/content/tools/intersect.html) | Find overlapping intervals in various ways. |
| [window](http://bedtools.readthedocs.org/en/latest/content/tools/window.html) | Find overlapping intervals within a window around an interval. |
| [closest](http://bedtools.readthedocs.org/en/latest/content/tools/closest.html) | Find the closest, potentially non-overlapping interval. |
| [coverage](http://bedtools.readthedocs.org/en/latest/content/tools/coverage.html) | Compute the coverage over defined intervals. |
| [map](http://bedtools.readthedocs.org/en/latest/content/tools/map.html) | Apply a function to a column for each overlapping interval. |
| [genomecov](http://bedtools.readthedocs.org/en/latest/content/tools/genomecov.html) | Compute the coverage over an entire genome. |
| [merge](http://bedtools.readthedocs.org/en/latest/content/tools/merge.html) | Combine overlapping/nearby intervals into a single interval. |
| [cluster](http://bedtools.readthedocs.org/en/latest/content/tools/cluster.html) | Cluster (but don't merge) overlapping/nearby intervals. |
| [complement](http://bedtools.readthedocs.org/en/latest/content/tools/complement.html) | Extract intervals \_not\_ represented by an interval file. |
| [shift](http://bedtools.readthedocs.org/en/latest/content/tools/shift.html) | Adjust the position of intervals. |
| [subtract](http://bedtools.readthedocs.org/en/latest/content/tools/subtract.html) | Remove intervals based on overlaps b/w two files. |
| [slop](http://bedtools.readthedocs.org/en/latest/content/tools/slop.html) | Adjust the size of intervals. |
| [flank](http://bedtools.readthedocs.org/en/latest/content/tools/flank.html) | Create new intervals from the flanks of existing intervals. |
| [sort](http://bedtools.readthedocs.org/en/latest/content/tools/sort.html) | Order the intervals in a file. |
| [random](http://bedtools.readthedocs.org/en/latest/content/tools/random.html) | Generate random intervals in a genome. |
| [shuffle](http://bedtools.readthedocs.org/en/latest/content/tools/shuffle.html) | Randomly redistribute intervals in a genome. |
| [sample](http://bedtools.readthedocs.org/en/latest/content/tools/sample.html) | Sample random records from file using reservoir sampling. |
| [spacing](http://bedtools.readthedocs.org/en/latest/content/tools/spacing.html) | Report the gap lengths between intervals in a file. |
| [annotate](http://bedtools.readthedocs.org/en/latest/content/tools/annotate.html) | Annotate coverage of features from multiple files. |
| **Multi-way** |
| [multiinter](http://bedtools.readthedocs.org/en/latest/content/tools/multiinter.html) | Identifies common intervals among multiple interval files. |
| [unionbedg](http://bedtools.readthedocs.org/en/latest/content/tools/unionbedg.html) | Combines coverage intervals from multiple BEDGRAPH files. |
| **Paired-end** |
| [pairtobed](http://bedtools.readthedocs.org/en/latest/content/tools/pairtobed.html) | Find pairs that overlap intervals in various ways. |
| [pairtopair](http://bedtools.readthedocs.org/en/latest/content/tools/pairtopair.html) | Find pairs that overlap other pairs in various ways. |
| **Format** |
| [bamtobed](http://bedtools.readthedocs.org/en/latest/content/tools/bamtobed.html) | Convert BAM alignments to BED (& other) formats. |
| [bedtobam](http://bedtools.readthedocs.org/en/latest/content/tools/bedtobam.html) | Convert intervals to BAM records. |
| [bamtofastq](http://bedtools.readthedocs.org/en/latest/content/tools/bamtofastq.html) | Convert BAM records to FASTQ records. |
| [bedpetobam](http://bedtools.readthedocs.org/en/latest/content/tools/bedpetobam.html) | Convert BEDPE intervals to BAM records. |
| [bed12tobed6](http://bedtools.readthedocs.org/en/latest/content/tools/bed12tobed6.html) | Breaks BED12 intervals into discrete BED6 intervals. |
| **Fasta** |
| [getfasta](http://bedtools.readthedocs.org/en/latest/content/tools/getfasta.html) | Use intervals to extract sequences from a FASTA file. |
| [maskfasta](http://bedtools.readthedocs.org/en/latest/content/tools/maskfasta.html) | Use intervals to mask sequences from a FASTA file. |
| [nuc](http://bedtools.readthedocs.org/en/latest/content/tools/nuc.html) | Profile the nucleotide content of intervals in a FASTA file. |
| **BAM** |
| [multicov](http://bedtools.readthedocs.org/en/latest/content/tools/multicov.html) | Counts coverage from multiple BAMs at specific intervals. |
| [tag](http://bedtools.readthedocs.org/en/latest/content/tools/tag.html) | Tag BAM alignments based on overlaps with interval files. |
| **Statistical** |
| [jaccard](http://bedtools.readthedocs.org/en/latest/content/tools/jaccard.html) | Calculate the Jaccard statistic b/w two sets of intervals. |
| [reldist](http://bedtools.readthedocs.org/en/latest/content/tools/reldist.html) | Calculate the distribution of relative distances b/w two files. |
| [fisher](http://bedtools.readthedocs.org/en/latest/content/tools/fisher.html) | Calculate Fisher statistic b/w two feature files. |
| **Miscellaneous** |
| [overlap](http://bedtools.readthedocs.org/en/latest/content/tools/overlap.html) | Computes the amount of overlap from two intervals. |
| [igv](http://bedtools.readthedocs.org/en/latest/content/tools/igv.html) | Create an IGV snapshot batch script. |
| [links](http://bedtools.readthedocs.org/en/latest/content/tools/links.html) | Create a HTML page of links to UCSC locations. |
| [makewindows](http://bedtools.readthedocs.org/en/latest/content/tools/makewindows.html) | Make interval "windows" across a genome. |
| [groupby](http://bedtools.readthedocs.org/en/latest/content/tools/groupby.html) | Group by common cols. & summarize oth. cols. (~ SQL "groupBy") |
| [expand](http://bedtools.readthedocs.org/en/latest/content/tools/expand.html) | Replicate lines based on lists of values in columns. |
| [split](http://bedtools.readthedocs.org/en/latest/content/tools/split.html) | Split a file into multiple files with equal records or base pairs. |
| [summary](http://bedtools.readthedocs.org/en/latest/content/tools/summary.html) | Statistical summary of intervals in a file. |



### References:


* Aaron R. Quinlan and Ira M. Hall. BEDTools: a flexible suite of utilities
 for comparing genomic features. Bioinformatics 2010, 26:841-842
 [Pubmed](http://www.ncbi.nlm.nih.gov/pubmed/20110278) | 
 [PMC](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2832824/) | 
 [Journal](http://bioinformatics.oxfordjournals.org/content/26/6/841.long)


Documentation
* [Official
 documentation](http://bedtools.readthedocs.org/en/latest/)
* [GitHub repository](https://github.com/arq5x/bedtools2)
* [Mailing
 list](http://groups.google.com/group/bedtools-discuss)
* [Tutorial](http://quinlanlab.org/tutorials/bedtools/bedtools.html)


Important Notes
* Module Name: bedtools (see [the modules page](/apps/modules.html) for more information)
* Example files in `$BEDTOOLS_TEST_DATA`



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

[user@cn3144 ~]$ **module load bedtools**
[user@cn3144 ~]$ **cp $BEDTOOLS\_TEST\_DATA/\* .**
[user@cn3144 ~]$ **ls -lh**
total 250M
-rw-r--r-- 1 user group  115 Feb  3 09:43 chroms.txt
-rw-r--r-- 1 user group 245M Feb  3 09:43 ENCFF001KOM.bam
-rw-r--r-- 1 user group 2.1M Feb  3 09:43 ENCFF001XYE.bed
-rw-r--r-- 1 user group 2.2M Feb  3 09:43 mm9_knownGene_promoters.bed

```

Merge identical intervals. Names (column 4 of bed file) of
merged intervals are concatenated. All intervals are 2000nts, each
centered on a promoter from the knownGenes table. Note that all
example files used here are sorted identically. Note the negative
distance with indicates the number of bp overlap required
for merging.



```

[user@cn3144 ~]$ **bedtools merge -d -2000 -c 4 -o distinct \
 -i mm9\_knownGene\_promoters.bed \
 > mm9\_knownGene\_upromoters.bed**
[user@cn3144 ~]$ **wc -l mm9\_knownGene\_promoters.bed**
55105 mm9_knownGene_promoters.bed
[user@cn3144 ~]$ **wc -l mm9\_knownGene\_upromoters.bed**
38862 mm9_knownGene_upromoters.bed

```

For each promoter interval, add a new column that contains the count of
H3K27ac peaks that overlap it

```

[user@cn3144 ~]$ **bedtools intersect -c -sorted \
 -a mm9\_knownGene\_upromoters.bed \
 -b ENCFF001XYE.bed \
 > mm9\_knownGene\_upromoters\_ac.bed**
[user@cn3144 ~]$ **head -5 mm9\_knownGene\_upromoters\_ac.bed**
chr1    3204713 3206713 uc007aet.1:     0
chr1    3647985 3649985 uc007aev.1:     0
chr1    3660579 3662579 uc007aeu.1:Q5GH67       0
chr1    4349395 4351395 uc007aex.2:Q548Q8       0
chr1    4398322 4400322 uc007aew.1:NP_001182591 0


```

Then count the number of reads from the H3K27ac IP that overlap the
same interval



```

[user@cn3144 ~]$ **bedtools intersect -c -sorted \
 -a mm9\_knownGene\_upromoters\_ac.bed \
 -b ENCFF001KOM.bam \
 > mm9\_knownGene\_upromoters\_ac\_reads.bed**

```

Show the mean and median number of reads per interval for promoters that overlap
0, 1, ... H3K27ac peaks. File needs to be sorted by the column that is used to
group lines



```

[user@cn3144 ~]$ **sort -k5,5 mm9\_knownGene\_upromoters\_ac\_reads.bed \
 | bedtools groupby -g 5 -c 6,6,6 -o count,mean,median**
0       24848   7.26295878943979378306  3
1       12288   95.410400390625 67
2       1675    89.6119402985074628987  74
3       48      86.4791666666666714036  62
4       3       43.6666666666666642982  49

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. bedtools.sh), which uses the input file 'bedtools.in'. For example:



```

#!/bin/bash 
set -e
set -o pipefail

module load bedtools
bedtools merge -d -2000 -c 4 -o distinct \
    -i mm9_knownGene_promoters.bed\
  > mm9_knownGene_upromoters.bed
bedtools intersect -c -sorted \
    -a mm9_knownGene_upromoters.bed \
    -b ENCFF001XYE.bed \
    > mm9_knownGene_upromoters_ac.bed
bedtools intersect -c -sorted \
    -a mm9_knownGene_upromoters_ac.bed \
    -b  ENCFF001KOM.bam \
    > mm9_knownGene_upromoters_ac_reads.bed
sort -k5,5 mm9_knownGene_upromoters_ac_reads.bed \
  | bedtools groupby -g 5 -c 6,6,6 -o count,mean,median \
  > median_counts.tsv

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] bedtools.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. bedtools.swarm). For example:



```

bamToBed -i file1.bam > file1.bed
bamToBed -i file2.bam > file2.bed
bamToBed -i file3.bam > file3.bed

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f bedtools.swarm [-g #] [-t #] --module bedtools
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g #  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t #  Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module bedtools  Loads the bedtools module for each subjob in the swarm 
 | |
 | |
 | |












