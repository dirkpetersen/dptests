

document.querySelector('title').textContent = 'trinity on Biowulf';
trinity on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
 |


 Trinity represents a novel method for the efficient and robust de novo
reconstruction of transcriptomes from RNA-seq data. Trinity combines three
independent software modules: Inchworm, Chrysalis, and Butterfly, applied
sequentially to process large volumes of RNA-seq reads. Trinity partitions the
sequence data into many individual de Bruijn graphs, each representing the
transcriptional complexity at at a given gene or locus, and then processes each
graph independently to extract full-length splicing isoforms and to tease apart
transcripts derived from paralogous genes. Briefly, the process works like so:


* Inchworm assembles the RNA-seq data into the unique sequences of
 transcripts, often generating full-length transcripts for a dominant
 isoform, but then reports just the unique portions of alternatively spliced
 transcripts.
* Chrysalis clusters the Inchworm contigs into clusters and constructs
 complete de Bruijn graphs for each cluster. Each cluster represents the
 full transcriptonal complexity for a given gene (or sets of genes that
 share sequences in common). Chrysalis then partitions the full read set
 among these disjoint graphs.
* Butterfly then processes the individual graphs in parallel, tracing the
 paths that reads and pairs of reads take within the graph, ultimately
 reporting full-length transcripts for alternatively spliced isoforms, and
 teasing apart transcripts that corresponds to paralogous genes.


 Trinity was developed at the Broad Institute & the Hebrew University of
Jerusalem.



In addition to these core functions, Trinity also incudes scripts to 
do in silico normalization, transcript quantitation, differential expression,
and other downstream analyses.


 Trinotate, the comprehensive annotation suite designed for automatic
functional annotation of transcriptomes, particularly de novo assembled
transcriptomes, from model or non-model organisms, is also available.


### References:


* Manfred G. Grabherr et al. *Trinity: reconstructing a full-length 
 transcriptome without a genome from RNA-Seq data*. Nature Biotechnology 
 2011, 29:644-652.
 [Pubmed](http://www.ncbi.nlm.nih.gov/pubmed/21572440) | 
 [PMC](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3571712/) | 
 [Journal](http://www.nature.com/nbt/journal/v29/n7/abs/nbt.1883.html)


Documentation
* [GitHub](https://github.com/trinityrnaseq/trinityrnaseq)
* [Manual](https://github.com/trinityrnaseq/trinityrnaseq/wiki)


Important Notes
* Module Name: trinity (see [the modules page](/apps/modules.html) for more information)
* Trinity is a multithreaded application. Cluster mode is not supported on biowulf any longer
* For some subset of trinity functions like, for example, expression estimation additional
 modules may have to be loaded.


Trinity creates a lot of temporary files and for efficiency we highly
recomment that it be run from **lscratch** as shown in the examples
below. This is especially true for swarms of trinity runs which can result in
severe stress on the shared file systems if lscratch is not used.


Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --cpus-per-task=6 --gres=lscratch:150 --mem=20g**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144]$ **module load trinity**
[user@cn3144]$ **cd /lscratch/$SLURM\_JOB\_ID**
[user@cn3144]$ **cp -pr $TRINITY\_TEST\_DATA/test\_Trinity\_Assembly .**
[user@cn3144]$ **cd test\_Trinity\_Assembly**
[user@cn3144]$ **Trinity --seqType fq --max\_memory 2G \
 --left reads.left.fq.gz \
 --right reads.right.fq.gz \
 --SS\_lib\_type RF \
 --CPU 4**
[...snip...]

[user@cn3144]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).

A Trinity run can be devided into two phases - (I) one phase that uses more memory
but is not easily parallelized and (II) a phase that uses less memory but can be
parallelized. Below is a trace of the memory usage and the number of active
threads for a genome based Trinity (2.4.0) assembly of ~15M normalized paired end human reads 
running on a single node with 10 allocated CPUs:




![Trinity memory/thread profile](/images/trinity_profile_singlenode.png)

While Trinity can parallelize Phase II across multiple nodes, we currently
do not support this capability on biowulf.


To run Trinity on a single node, create a batch script similar to the
following example.



```

#! /bin/bash
# this file is trinity.sh
function die() {
    echo "$@" >&2
    exit 1
}
module load trinity/2.14.0 || die "Could not load trinity module"
[[ -d /lscratch/$SLURM_JOB_ID ]] || die "no lscratch allocated"

inbam=$1
mkdir /lscratch/$SLURM_JOB_ID/in
mkdir /lscratch/$SLURM_JOB_ID/out
cp $inbam /lscratch/$SLURM_JOB_ID/in
bam=/lscratch/$SLURM_JOB_ID/in/$(basename $inbam)
out=/lscratch/$SLURM_JOB_ID/out

Trinity --genome_guided_bam $bam \
    --SS_lib_type RF \
    --output  $out \
    --genome_guided_max_intron 10000 \
    --max_memory 28G \
    --CPU 12
mv $out/Trinity-GG.fasta /data/$USER/trinity_out/$(basename $inbam .bam)-Trinity-GG.fasta

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```

biowulf$ **sbatch --mem=30g --cpus-per-task=12 --gres=lscratch:150 trinity.sh /data/$USER/trinity\_in/sample.bam**

```





