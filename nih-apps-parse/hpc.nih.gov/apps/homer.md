

document.querySelector('title').textContent = 'homer on Biowulf';
homer on Biowulf


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


homer solid tools and methods for interpreting ChIP-Seq and other HTS
experiments. In addition to UCSC visualization support, peak finding, and
motif finding of course, homer can help assemble data across multiple
experiments and look at position-specific relationships between sequencing
tags, motifs, and other features. You do not need to use the peak finding
methods in this package to use motif finding. (Use the bed2pos.pl program to
create peak files from BED files).


### References:


* S. Heinz, C. Benner, N. Spann, E. Bertolino et al. 
 *Simple Combinations of Lineage-Determining Transcription Factors Prime cis-Regulatory Elements
 Required for Macrophage and B Cell Identities*. Mol Cell 2010(4):576-589.
 [PubMed](https://www.ncbi.nlm.nih.gov/pubmed/20513432) | 
 [PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2898526/) | 
 [Journal](https://www.sciencedirect.com/science/article/pii/S1097276510003667?via%3Dihub)


Documentation
* [Home](http://homer.ucsd.edu/homer)


Important Notes
* Module Name: homer (see [the modules page](/apps/modules.html) for more information)
* Example files in `$HOMER_TEST_DATA`
* Reference data in /fdb/homer/



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=10g --gres=lscratch:20**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load homer**
[user@cn3144 ~]$ **cd /lscratch/$SLURM\_JOB\_ID**
[user@cn3144 ~]$ # copy some test data (CTCF ChIP-Seq from human Parathyroid adenoma)
[user@cn3144 ~]$ **cp -L $HOMER\_TEST\_DATA/\* .**
[user@cn3144 ~]$ **ls -lh**
-rw-rw-r-- 1 user group 869M Mar 15 10:42 ENCFF374UMR.bam
[user@cn3144 ~]$ **makeTagDirectory ctcf -unique -genome hg38 -checkGC ENCFF374UMR.bam**
        Will parse file: ENCFF374UMR.bam

        Creating directory: ctcf and removing existing *.tags.tsv

        Treating ENCFF374UMR.bam as a bam file
        Reading alignment file ENCFF374UMR.bam

        Optimizing tag files...
        Estimated genome size = 3098172877
        Estimated average read density = 0.017473 per bp
        Total Tags = 54133547.0
        Total Positions = 43151510
        Average tag length = 76.0
        Median tags per position = 1 (ideal: 1)
        Average tags per position = 1.254
        Fragment Length Estimate: 198
        Peak Width Estimate: 333
        Autocorrelation quality control metrics:
                Same strand fold enrichment: 2.0
                Diff strand fold enrichment: 1.9
                Same / Diff fold enrichment: 1.0
        [...snip...]
[user@cn3144 ~]$ **du -sh ctcf**
1021M   ctcf
[user@cn3144 ~]$ **ml ucsc**
[user@cn3144 ~]$ **samtools view -H ENCFF374UMR.bam | grep '@SQ' | tr ':' '\t' | cut -f3,5 > chroms**
[user@cn3144 ~]$ **makeUCSCfile ctcf -bigWig chroms -fsize 1e20 -o ctcf.bw**
...

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. homer.sh), which uses the input file 'homer.in'. For example:



```

#!/bin/bash
module load homer/4.9.1 || exit 1
findPeaks ctcf_tagdir -style factor -o auto -i control_tags

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=2 --mem=10g homer.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. homer.swarm). For example:



```

cd /data/$USER/dir1; findMotifs.pl genelist.txt mouse Results -len 10
cd /data/$USER/dir2; findMotifs.pl genelist.txt mouse Results -len 10
cd /data/$USER/dir3; findMotifs.pl genelist.txt mouse Results -len 10

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f homer.swarm -g 10 -t 2 --module homer/4.9.1
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g #  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t #  Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module homer  Loads the homer module for each subjob in the swarm 
 | |
 | |
 | |








