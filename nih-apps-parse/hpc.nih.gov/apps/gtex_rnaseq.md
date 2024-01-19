

document.querySelector('title').textContent = 'gtex\_rnaseq on Biowulf';
gtex\_rnaseq on Biowulf


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


 This module makes available the tools used in the [GTEX
 RNA-Seq pipeline](https://github.com/broadinstitute/gtex-pipeline/tree/master/rnaseq). Also planned is the implementation of a combined
 pipeline but that is not yet ready.



Documentation
* gtex\_rnaseq on [GitHub](https://github.com/broadinstitute/gtex-pipeline/tree/master/rnaseq)


Important Notes
* Module Name: gtex\_rnaseq (see [the modules page](/apps/modules.html) for more information)
* Some of the tools are multithreaded. When running individual tools please match the number
 of allocated CPUs to the number of threads
* Reference data in `GTEX_RNASEQ_REF`
* Example files in `GTEX_RNASEQ_TEST_DATA`



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the
individual steps of the pipeline. Note that only some of the steps are resource
intensive. For production work please run individual steps as batch jobs with
the appropriate resource allocations. This full size example takes some hours
to run even though it is downsampled to just 10M reads per sample.


Note that in this interactive examples all steps are run as part of the same allocation.
That is ok for testing and debugging, but since the steps require quite different resources
(e.g. star scales well to multiple threads and requires around 45GB of memory, rnaseqc and
markduplicates are mostly single threaded) for production each step should be run as a separate
job.



```

[user@biowulf]$ **sinteractive --mem=40g --time=1-12 --cpus-per-task=16 --gres=lscratch:400**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144]$ **module load gtex\_rnaseq**
# temp files for STAR are located in the same directory as output so better be in lscratch
[user@cn3144]$ **cd /lscratch/$SLURM\_JOB\_ID**
[user@cn3144]$ **cp -Lr ${GTEX\_RNASEQ\_TEST\_DATA:-none} ./data**
[user@cn3144]$ **samples=( HepG2\_ENCLB059ZZZ 
 HepG2\_ENCLB060ZZZ 
 K562\_ENCLB063ZZZ 
 K562\_ENCLB064ZZZ )**
[user@cn3144]$  **for sample in "${samples[@]}" ; do
 run\_STAR.py \
 --output\_dir data/star\_out \
 -t $SLURM\_CPUS\_PER\_TASK \
 ${GTEX\_RNASEQ\_REF}/GRCh38\_release38/star\_index\_oh100 \
 data/${sample}\_R1.fastq.gz \
 data/${sample}\_R2.fastq.gz \
 ${sample}
 done**

# note that this script changes into the output directory,
# so the path to the input bam has to be absolute or relative to
# that output directory
[user@cn3144]$ **for sample in "${samples[@]}" ; do
 run\_MarkDuplicates.py \
 -o data \
 star\_out/${sample}.Aligned.sortedByCoord.out.bam \
 ${sample}.Aligned.sortedByCoord.out.md
 done**

# need to manually index the markduplicate bam output
[user@cn3144]$ **module load samtools**
[user@cn3144]$ **for sample in "${samples[@]}" ; do
 samtools index data/${sample}.Aligned.sortedByCoord.out.md.bam
 done**

# - Syntax for the V8 pipeline:
#   - need to use the older java
#   - need more memory
[user@cn3144]$ **for sample in "${samples[@]}" ; do
 run\_rnaseqc.py \
 --output\_dir=$PWD/data \
 --java /usr/lib/jvm/java-1.7.0-openjdk-amd64/bin/java \
 --memory=16 \
 ${sample}.Aligned.sortedByCoord.out.md.bam \
 ${GTEX\_RNASEQ\_REF}/GRCh38\_release38/gencode.v38.primary\_assembly.genes.gtf \
 ${GTEX\_RNASEQ\_REF}/GRCh38\_release38/GRCh38.primary\_assembly.genome.fa \
 ${sample}
 done**

# - Syntax for the V10 pipeline. Need to change to the directory for this
#   step to work
[user@cn3144]$ **pushd data**
[user@cn3144]$ **for sample in "${samples[@]}" ; do
 run\_rnaseqc.py \
 ${GTEX\_RNASEQ\_REF}/GRCh38\_release38/gencode.v38.primary\_assembly.genes.gtf \
 ${sample}.Aligned.sortedByCoord.out.md.bam \
 ${sample}
 done**
[user@cn3144]$ **popd** 

[user@cn3144]$ **for sample in "${samples[@]}" ; do
 run\_RSEM.py \
 --threads $SLURM\_CPUS\_PER\_TASK \
 ${GTEX\_RNASEQ\_REF}/GRCh38\_release38/rsem\_reference \
 $PWD/data/star\_out/${sample}.Aligned.toTranscriptome.out.bam \
 $PWD/data/${sample}
 done**

[user@cn3144]$ **combine\_GCTs.py data/\*.exon\_reads.gct.gz data/combined\_gcts\_exon\_reads**

# copy relevant results back

[user@cn3144]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. gtex\_rnaseq.sh), which uses the input file 'gtex\_rnaseq.in'. For example:



```

#!/bin/bash
module load gtex_rnaseq/V8 || exit 1
wd=$PWD
sample=HepG2_ENCLB059ZZZ
cd /lscratch/$SLURM_JOB_ID || exit 1
mkdir data
cp -L ${GTEX_RNASEQ_TEST_DATA:-none}/${sample}* ./data
run_STAR.py \
    --output_dir $PWD/data/star_out \
    -t $SLURM_CPUS_PER_TASK \
    ${GTEX_RNASEQ_REF}/GRCh38_release38/star_index_oh100 \
    data/${sample}_R1.fastq.gz \
    data/${sample}_R2.fastq.gz \
    ${sample}
cp -r data $wd

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=16 --mem=40g --gres=lscratch:75 gtex_rnaseq.sh
```







