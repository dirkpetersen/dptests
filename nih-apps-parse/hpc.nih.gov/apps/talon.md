

document.querySelector('title').textContent = 'talon on Biowulf';
talon on Biowulf


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




>  TALON is a Python package for identifying and quantifying known
>  and novel genes/isoforms in long-read transcriptome data sets. TALON is
>  technology-agnostic in that it works from mapped SAM files, allowing data
>  from different sequencing platforms (i.e. PacBio and Oxford Nanopore) to be
>  analyzed side by side. 


### References:


* Dana Wyman et al., *A technology-agnostic long-read analysis pipeline for transcriptome discovery and quantification*. [bioRxiv](https://doi.org/10.1101/672931)



Documentation
* talon on [GitHub](https://github.com/mortazavilab/TALON)


Important Notes
* Module Name: talon (see [the modules page](/apps/modules.html) for more information)
* Some talon tools are multithreaded. Please match the number of threads with the number
 of allocated CPUs
* Example files in `$TALON_TEST_DATA`
* Some talon steps may require large amounts of memory. See [this GitHub issue](https://github.com/mortazavilab/TALON/issues/90) for a discussion on how to reduce memory usage.



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run through
the steps using 2 replicates of human cardiac atrium tissue runs on a PacBio Sequel II:



```

[user@biowulf]$ **sinteractive --cpus-per-task=6 --mem=16G --gres=lscratch:50**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144]$ **cd /lscratch/$SLURM\_JOB\_ID**
[user@cn3144]$ **module load talon**
[user@cn3144]$ **cp -L ${TALON\_TEST\_DATA:-none}/\* .**
[user@cn3144]$ **ls -lh**
total 4.6G
-rw-r--r-- 1 user group 172M Apr 20 14:31 ENCFF291EKY.bam
-rw-r--r-- 1 user group 1.7M Apr 20 14:31 ENCFF291EKY.bam.bai
-rw-r--r-- 1 user group 189M Apr 20 14:31 ENCFF613SDS.bam
-rw-r--r-- 1 user group 1.7M Apr 20 14:31 ENCFF613SDS.bam.bai
-rw-r--r-- 1 user group 1.3G Apr 20 14:31 gencode.v35.primary_assembly.annotation.gtf
-rw-r--r-- 1 user group 3.0G Apr 20 14:31 GRCh38.primary_assembly.genome.fa
-rw-r--r-- 1 user group 6.4K Apr 20 14:31 GRCh38.primary_assembly.genome.fa.fai
[user@cn3144]$ **gtf=gencode.v35.primary\_assembly.annotation.gtf**
[user@cn3144]$ **genome=GRCh38.primary\_assembly.genome.fa**
[user@cn3144]$ **bam1=ENCFF291EKY.bam**
[user@cn3144]$ **bam2=ENCFF613SDS.bam**
[user@cn3144]$ **talon\_initialize\_database \
 --f $gtf \
 --a gencode\_35 \
 --g GRCh38 \
 --o example\_talon**
chr1
bulk update genes...
bulk update gene_annotations...
bulk update transcripts...
[...snip...]
[user@cn3144]$ **mkdir -p labeled tmp**
[user@cn3144]$ ### check internal priming sites
[user@cn3144]$ **talon\_label\_reads --f $bam1 \
 --g $genome \
 --t $SLURM\_CPUS\_PER\_TASK \
 --ar 20 \
 --tmpDir=/lscratch/$SLURM\_JOB\_ID/tmp \
 --deleteTmp \
 --o labeled/${bam1%.bam}**
[ 2021-04-20 17:10:44 ] Started talon_label_reads run.
[ 2021-04-20 17:10:44 ] Splitting SAM by chromosome...
[ 2021-04-20 17:10:44 ] -----Writing chrom files...
[ 2021-04-20 17:10:59 ] Launching parallel jobs...
[ 2021-04-20 17:11:14 ] Pooling output files...
[ 2021-04-20 17:11:27 ] Run complete
[user@cn3144]$ **talon\_label\_reads --f $bam2 \
 --g $genome \
 --t $SLURM\_CPUS\_PER\_TASK \
 --ar 20 \
 --tmpDir=/lscratch/$SLURM\_JOB\_ID/tmp \
 --deleteTmp \
 --o labeled/${bam2%.bam}**
[...snip...]
[user@cn3144]$ ### run talon annotator
[user@cn3144]$ **cat > config.csv <<\_\_EOF\_\_**
ex_rep1,GRCh38,PacBio-Sequel2,labeled/${bam1%.bam}_labeled.sam
ex_rep2,GRCh38,PacBio-Sequel2,labeled/${bam2%.bam}_labeled.sam
**\_\_EOF\_\_**
[user@cn3144]$ **talon \
 -t $SLURM\_CPUS\_PER\_TASK \
 --f config.csv \
 --db example\_talon.db \
 --build GRCh38 \
 --o example**
[user@cn3144]$ ### summarize results
[user@cn3144]$  **talon\_summarize \
 --db example\_talon.db \
 --v \
 --o example**
[user@cn3144]$ ### run any other tools and then copy results back to shared space
[user@cn3144]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. talon.sh), which uses the input file 'talon.in'. For example:



```

#! /bin/bash

module load talon/5.0

bam1=ENCFF291EKY.bam
bam2=ENCFF613SDS.bam
gtf=gencode.v35.primary_assembly.annotation.gtf
genome=GRCh38.primary_assembly.genome.fa



cd /lscratch/$SLURM_JOB_ID
cp -L ${TALON_TEST_DATA:-none}/* .
talon_initialize_database \
    --f $gtf \
    --a gencode_35 \
    --g GRCh38 \
    --o example_talon

mkdir -p labeled tmp
talon_label_reads --f $bam1 \
    --g $genome  \
    --t $SLURM_CPUS_PER_TASK \
    --ar 20 \
    --tmpDir=/lscratch/$SLURM_JOB_ID/tmp \
    --deleteTmp \
    --o labeled/${bam1%.bam}
talon_label_reads --f $bam2 \
    --g $genome  \
    --t $SLURM_CPUS_PER_TASK \
    --ar 20 \
    --tmpDir=/lscratch/$SLURM_JOB_ID/tmp \
    --deleteTmp \
    --o labeled/${bam2%.bam}

cat > config.csv <<__EOF__
ex_rep1,GRCh38,PacBio-Sequel2,labeled/${bam1%.bam}_labeled.sam
ex_rep2,GRCh38,PacBio-Sequel2,labeled/${bam2%.bam}_labeled.sam
__EOF__

talon \
    -t $SLURM_CPUS_PER_TASK \
    --f config.csv \
    --db example_talon.db \
    --build GRCh38 \
    --o example

talon_summarize \
   --db example_talon.db \
   --v \
   --o example

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] talon.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. talon.swarm). For example:



```

talon_label_reads --f ENCFF291EKY.bam \
    --g GRCh38.primary_assembly.genome.fa  \
    --t $SLURM_CPUS_PER_TASK \
    --ar 20 \
    --tmpDir=/lscratch/$SLURM_JOB_ID \
    --deleteTmp \
    --o labeled/ENCFF291EKY
talon_label_reads --f ENCFF613SDS.bam \
    --g GRCh38.primary_assembly.genome.fa  \
    --t $SLURM_CPUS_PER_TASK \
    --ar 20 \
    --tmpDir=/lscratch/$SLURM_JOB_ID \
    --deleteTmp \
    --o labeled/ENCFF613SDS

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f talon.swarm [-g 10] [-t 6] --gres=lscratch:50 --module talon/5.0
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g #  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t #  Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module talon  Loads the talon module for each subjob in the swarm 
 | |
 | |
 | |








