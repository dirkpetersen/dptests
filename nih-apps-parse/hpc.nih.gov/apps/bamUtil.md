

document.querySelector('title').textContent = 'Bamutil on Biowulf';
Bamutil on Biowulf


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



bamUtil is a collection of utilities for manipulating bam files compiled
into a single executable called `bam`. Available tools:





| **Tools to rewrite SAM/BAM files** |
| --- |
|  [convert](http://genome.sph.umich.edu/wiki/BamUtil:_convert )  |  Convert SAM/BAM to SAM/BAM |
|  [writeRegion](http://genome.sph.umich.edu/wiki/BamUtil:_writeRegion )  |  Write a file with reads in the specified region and/or have the specified read name |
|  [splitChromosome](http://genome.sph.umich.edu/wiki/BamUtil:_splitChromosome )  |  Split BAM by Chromosome |
|  [splitBam](http://genome.sph.umich.edu/wiki/BamUtil:_splitBam )  |  Split a BAM file into multiple BAM files based on ReadGroup |
|  [findCigars](http://genome.sph.umich.edu/wiki/BamUtil:_findCigars )  |  Output just the reads that contain any of the specified CIGAR operations. |
| **Tools to modify SAM/BAM files** |
|  [clipOverlap](http://genome.sph.umich.edu/wiki/BamUtil:_clipOverlap )  |  Clip overlapping read pairs in a SAM/BAM File already sorted by Coordinate or ReadName |
|  [filter](http://genome.sph.umich.edu/wiki/BamUtil:_filter )  |  Filter reads by clipping ends with too high of a mismatch percentage and by marking reads unmapped if the quality of mismatches is too high |
|  [revert](http://genome.sph.umich.edu/wiki/BamUtil:_revert )  |  Revert SAM/BAM replacing the specified fields with their previous values (if known) and removes specified tags |
|  [squeeze](http://genome.sph.umich.edu/wiki/BamUtil:_squeeze )  |  reduces files size by dropping OQ fields, duplicates, & specified tags, using '=' when a base matches the reference, binning quality scores, and replacing readNames with unique integers |
|  [trimBam](http://genome.sph.umich.edu/wiki/BamUtil:_trimBam )  |  Trim the ends of reads in a SAM/BAM file changing read ends to 'N' and quality to '!' |
|  [mergeBam](http://genome.sph.umich.edu/wiki/BamUtil:_mergeBam )  |  merge multiple BAMs and headers appending ReadGroupIDs if necessary |
|  [polishBam](http://genome.sph.umich.edu/wiki/BamUtil:_polishBam )  |  adds/updates header lines & adds the RG tag to each record |
|  [dedup](http://genome.sph.umich.edu/wiki/BamUtil:_dedup )  |  Mark Duplicates |
|  [recab](http://genome.sph.umich.edu/wiki/BamUtil:_recab )  |  Recalibrate |
| **Tools to extract information** |
|  [validate](http://genome.sph.umich.edu/wiki/BamUtil:_validate )  |  Validate a SAM/BAM File |
|  [diff](http://genome.sph.umich.edu/wiki/BamUtil:_diff )  |  Diff 2 coordinate sorted SAM/BAM files. |
|  [stats](http://genome.sph.umich.edu/wiki/BamUtil:_stats )  |  Stats a SAM/BAM File |
|  [gapInfo](http://genome.sph.umich.edu/wiki/BamUtil:_gapInfo )  |  Print information on the gap between read pairs in a SAM/BAM File. |
|  [dumpHeader](http://genome.sph.umich.edu/wiki/BamUtil:_dumpHeader )  |  Print SAM/BAM Header |
|  [dumpRefInfo](http://genome.sph.umich.edu/wiki/BamUtil:_dumpRefInfo )  |  Print SAM/BAM Reference Name Information |
|  [dumpIndex](http://genome.sph.umich.edu/wiki/BamUtil:_dumpIndex )  |  Print BAM Index File in English |
|  [readReference](http://genome.sph.umich.edu/wiki/BamUtil:_readReference )  |  Print the reference string for the specified region |
|  [explainFlags](http://genome.sph.umich.edu/wiki/BamUtil:_explainFlags )  |  Describe flags |
|  [bam2FastQ](http://genome.sph.umich.edu/wiki/BamUtil:_bam2FastQ )  |  Convert the specified BAM file to fastQs. |
|  [readIndexedBam](http://genome.sph.umich.edu/wiki/BamUtil:_readIndexedBam )  |  Read Indexed BAM By Reference and write it from reference id  |


Documentation
* [bamutil Main Site](https://genome.sph.umich.edu/wiki/BamUtil#bamUtil_Overview)


Important Notes
* Module Name: bamutil (see [the modules page](/apps/modules.html) for more information)
* Singlethreaded app



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

[user@cn3144 ~]$ **module load bamutil**
[+] Loading bamutil  1.0.14 

[user@cn3144 ~]$ **cd /data/$USER/test\_data/bam**

[user@cn3144 ~]$ **bam stats --in read1\_250k.bam --qual --basic**
Number of records read = 170016

TotalReads      170016.00
MappedReads     170016.00
PairedReads     0.00
ProperPair      0.00
.....

[user@cn3144 ~]$ **bam explainFlags --dec 96**
0x60 (96):
        mate reverse strand
        first fragment

[user@cn3144 ~]$ **bam dedup --in read1\_500k\_sorted.bam --out /scratch/$USER/temp.bam \**
**--rmDups**

[user@cn3144 ~]$ **cat /scratch/$USER/temp.bam.log**
......
--------------------------------------------------------------------------
SUMMARY STATISTICS OF THE READS
Total number of reads: 326389
Total number of paired-end reads: 0
Total number of properly paired reads: 0
Total number of unmapped reads: 0
Total number of reverse strand mapped reads: 162220
Total number of QC-failed reads: 0
Total number of secondary reads: 0
Size of singleKeyMap (must be zero): 0
Size of pairedKeyMap (must be zero): 0
Total number of missing mates: 0
Total number of reads excluded from duplicate checking: 0
--------------------------------------------------------------------------
Sorting the indices of 3651 duplicated records

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. bamutil.sh). For example:



```

#!/bin/bash
#SBATCH --job-name=bamUtil

set -e
module load bamutil
inbam=/data/$USER/test_data/bam/gcat_set_025_raw.bam
obam=/data/$USER/test_data/temp/gcat_set_025.clean.bam
dbsnp=/data/$USER/test_data/snp/snp138_pos.hg19
genome=/fdb/igenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa

# bamutil wants to write to the directory containing the genome file
mkdir /lscratch/$SLURM_JOB_ID
cp $genome /lscratch/$SLURM_JOB_ID
lgenome=/lscratch/$SLURM_JOB_ID/$(basename $genome)

bam dedup --in $inbam --out $obam \
  --rmDups --verbose --recab \
  --refFile $lgenome \
  --dbsnp $dbsnp \
  --storeQualTag OQ --maxBaseQual 40
rm -rf /lscratch/$SLURM_JOB_ID

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --mem=4g bamutil.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. bamutil.swarm). For example:



```

bam squeeze --in file1.bam --out file1s.bam --refFile genome.fa 
bam squeeze --in file2.bam --out file2s.bam --refFile genome.fa 
bam squeeze --in file3.bam --out file3s.bam --refFile genome.fa 

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f bamutil.swarm --module bamutil
```

where


|  |  |
| --- | --- |
| --module bamutil Loads the bamutil module for each subjob in the swarm 
 | |








