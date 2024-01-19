

document.querySelector('title').textContent = 'STAR and STAR-Fusion on Biowulf';
STAR and STAR-Fusion on Biowulf


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


STAR aligns RNA-Seq data to reference genomes. It is designed to be fast and
accurate for known and novel splice junctions. In addition, it has no limit on
the read size and can align reads with multiple splice junctions. This is
becoming more important as read lengths increase. As a result of the alignment
algorithm, poor quality tails and poly-A tails are clipped in the alignment.


Starting with version 2.4.1a, annotations can be included at the mapping stage
instead/in addition to the indexing stage.


STAR-Fusion uses the STAR aligner to identify candidate fusion transcripts.


There are now separate modules for STAR and STAR-Fusion.


### References:


* Alexander Dobin *et al.*
*STAR: ultrafast universal RNA-seq aligner*.
 Bioinformatics 2013, 29:15-21.
 [Pubmed](http://www.ncbi.nlm.nih.gov/pubmed/23104886) | 
 [PMC](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3530905/) | 
 [Journal](http://bioinformatics.oxfordjournals.org/content/29/1/15.long)


Documentation
* [STAR Home page](https://github.com/alexdobin/STAR)
* [STAR Manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)
* [STAR Google group](https://groups.google.com/forum/#!forum/rna-star)
* [STAR-Fusion wiki](https://github.com/STAR-Fusion/STAR-Fusion/wiki)


Important Notes
* Module Name: STAR (see [the modules page](/apps/modules.html) for more information)
* STAR is a multithreaded application. Please match the number of allocated CPUs with the
number of threads
* Example files in `$STAR_TEST_DATA`



STAR Indices
STAR made a backwards incompatible change to it's index structure 
with version 2.4.2.


STAR indices can be found in


  **`/fdb/STAR_indices/[STAR VERSION]`**


For historical reasons there are also three symlinks pointing to 
the versioned subdirectories in `/fdb/STAR_indices`:


* `**/fdb/STAR** → /fdb/STAR_indices/2.4.0d`
* `**/fdb/STAR\_2.4.2a** → /fdb/STAR_indices/2.4.2a`
* `**/fdb/STAR\_current** → /fdb/STAR_indices/[newest STAR]`


STAR-Fusion
Various CTAT libraries are available at **`/fdb/CTAT`**. Some of
the libraries are directly at the top level. Starting from STAR-Fusion 1.7, CTAT libraries
in the correct format are available in version specific subdirectories such as



```

/fdb/CTAT/__genome_libs_StarFv1.7
/fdb/CTAT/__genome_libs_StarFv1.9
/fdb/CTAT/__genome_libs_StarFv1.10

```

Note that 1.11.0 and 1.12.0 are still compatible with \_\_genome\_libs\_StarFv1.10


The following table shows the minimum version of STAR for recent versions of STAR-Fusion:




| **STAR-Fusion** | **Minumum STAR version** |
| --- | --- |
| 1.12.0 | 2.7.8a |
| 1.11.0 | 2.7.8a |
| 1.10.0 | 2.7.8a |
| 1.9.1 | 2.7.2b |
| 1.7.0 | 2.7.2a |
| 1.6.0 | 2.7.0f |
| 1.5.0 | 2.6.1 |


Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
STAR is used to create genome indices as well as to align and map short reads
to the indexed genome. A GTF format annotation of transcripts can be provided
during indexing or, since version 2.4.1a, on the fly during mapping.


A simple example of indexing a small genome with annotation. If annotation is
provided, a overhang depending on the readlength to be used has to be provided
as well. In this example we use a small genome, so 30g is more than enough
memory. For example, for 100nt reads:


Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --cpus-per-task=12 --mem=30g --gres=lscratch:20**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load STAR**
[user@cn3144 ~]$ **mkdir -p indices/star100-EF4**
[user@cn3144 ~]$ **GENOME=/fdb/igenomes/Saccharomyces\_cerevisiae/Ensembl/EF4**
[user@cn3144 ~]$ **STAR \
 --runThreadN 2 \
 --runMode genomeGenerate \
 --genomeDir indices/star100-EF4 \
 --genomeFastaFiles $GENOME/Sequence/WholeGenomeFasta/genome.fa \
 --sjdbGTFfile $GENOME/Annotation/Genes/genes.gtf \
 --sjdbOverhang 99 \
 --genomeSAindexNbases 10**

```

However, a default value of 100 will work well in many cases. Note that
for the output of sorted bam files, the temp directory should be placed in lscratch.


Align single end 50nt RNA-Seq data to the mouse genome:



```

[user@cn3144 ~]$ **mkdir -p bam/rnaseq\_STAR**
[user@cn3144 ~]$ **GENOME=/fdb/STAR\_current/UCSC/mm10/genes-50**
[user@cn3144 ~]$ **STAR \
 --runThreadN 12 \
 --genomeDir $GENOME \
 --sjdbOverhang 50 \
 --readFilesIn $STAR\_TEST\_DATA/ENCFF138LJO\_1M.fastq.gz \
 --readFilesCommand zcat \
 --outSAMtype BAM SortedByCoordinate \
 --outTmpDir=/lscratch/$SLURM\_JOB\_ID/STARtmp \
 --outFileNamePrefix bam/rnaseq\_STAR/test**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```

To use STAR-Fusion, please load the separately available STAR-Fusion module



```

[user@3144 ~]$ **module load STAR-Fusion**

```

For more details on STAR-Fusion see the [STAR-Fusion wiki](https://github.com/STAR-Fusion/STAR-Fusion/wiki).



Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Batch scripts should use the `$SLURM_CPUS_PER_TASK` environment
variable to determine the number of threads to use. This variable is set by
slurm according to the number of requested cores. An example script for
aligning single end RNA-Seq data with STAR might look like the following:



```

#! /bin/bash
# this file is STAR.sh
set -o pipefail
set -e

function fail {
    echo "$@" >&2
    exit 1
}

module load samtools/1.6 || fail "could not load samtools module"
module load STAR         || fail "could not load STAR module"
cd /data/$USER/test_data || fail "no such directory"
mkdir -p bam/rnaseq_STAR
GENOME=/fdb/STAR_current/UCSC/mm10/genes-50
STAR \
    --runThreadN $SLURM_CPUS_PER_TASK \
    --genomeDir $GENOME \
    --sjdbOverhang 50 \
    --readFilesIn fastq/rnaseq_500k.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --outTmpDir=/lscratch/$SLURM_JOB_ID/STARtmp \
    --outFileNamePrefix bam/rnaseq_STAR/test

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command. STAR requires
30-45g for mammalian genomes. Human genomes generally require about 45GB.



```
sbatch --cpus-per-task=12 --mem=38g --gres=lscratch:20 STAR.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. STAR.swarm). For example:



```

cd /data/$USER/test_data \
  && mkdir -p bam/rnaseq_STAR \
  && STAR \
    --runThreadN $SLURM_CPUS_PER_TASK \
    --genomeDir /fdb/STAR_current/UCSC/mm10/genes-50 \
    --sjdbOverhang 50 \
    --readFilesIn fastq/rnaseq_500k.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --outTmpDir=/lscratch/$SLURM_JOB_ID/STARtmp \
    --outFileNamePrefix bam/rnaseq_STAR/test
cd /data/$USER/test_data \
  && mkdir -p bam/rnaseq_STAR2 \
  && STAR \
    --runThreadN $SLURM_CPUS_PER_TASK \
    --genomeDir /fdb/STAR_current/UCSC/mm10/genes-50 \
    --sjdbOverhang 50 \
    --readFilesIn fastq/rnaseq_250k.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --outTmpDir=/lscratch/$SLURM_JOB_ID/STARtmp \
    --outFileNamePrefix bam/rnaseq_STAR2/test

```

Note the use of line continuation. Submit this job using the 
[swarm](/apps/swarm.html) command.



```
swarm -f STAR.swarm -t 12 -g 38 --module STAR
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g #  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t #  Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module STAR  Loads the STAR module for each subjob in the swarm 
 | |
 | |
 | |








