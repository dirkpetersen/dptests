

document.querySelector('title').textContent = 'Hicexplorer on Biowulf';
Hicexplorer on Biowulf


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


Hicexplorer is a set of tools to process, normalize and visualize Hi-C data.



### References:


* Wolff J, Rabbani L, Gilsbach R, Richard G, Manke T, Backofen R, Grüning BA.
*[Galaxy HiCExplorer 3: a web server for reproducible Hi-C, capture Hi-C and single-cell Hi-C data analysis, quality control and visualization.](https://pubmed.ncbi.nlm.nih.gov/32301980/)* Nucleic Acids Res. 2020 Jul 2;48(W1):W177-W184.


Documentation
* Hicexplorer Main Site: [GitHub](https://github.com/deeptools/HiCExplorer)


Important Notes
* Module Name: hicexplorer (see [the modules page](/apps/modules.html) for more information)



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=8g --cpus-per-task=8 --gres=lscratch:100**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn4224 are ready for job

[user@cn4224 ~]$ **module load hicexplorer bowtie samtools**
[+] Loading hicexplorer  3.5.1  on cn4224 
[+] Loading singularity  3.7.0  on cn4224
[+] Loading bowtie  2-2.4.2 
[+] Loading samtools 1.11  ... 

[user@cn4224 ~]$ **cd /lscratch/${SLURM\_JOB\_ID}**

[user@cn4224 ~]$ **wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR395/005/SRR3950565/SRR3950565\_1.fastq.gz**
[user@cn4224 ~]$ **wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR395/005/SRR3950565/SRR3950565\_2.fastq.gz**
[user@cn4224 ~]$ **wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR395/006/SRR3950566/SRR3950566\_1.fastq.gz**
[user@cn4224 ~]$ **wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR395/006/SRR3950566/SRR3950566\_2.fastq.gz**
[user@cn4224 ~]$ **wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR395/009/SRR3950559/SRR3950559\_1.fastq.gz**
[user@cn4224 ~]$ **wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR395/009/SRR3950559/SRR3950559\_2.fastq.gz**
[user@cn4224 ~]$ **wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR395/000/SRR3950560/SRR3950560\_1.fastq.gz**
[user@cn4224 ~]$ **wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR395/000/SRR3950560/SRR3950560\_2.fastq.gz**

[user@cn4224 ~]$ **bowtie2 -x mm9\_index --threads 8 -U SRR3950565\_1.fastq.gz --reorder | samtools view -Shb - > SRR3950565\_1.bam**
[user@cn4224 ~]$ **bowtie2 -x mm9\_index --threads 8 -U SRR3950565\_2.fastq.gz --reorder | samtools view -Shb - > SRR3950565\_2.bam**
[user@cn4224 ~]$ **bowtie2 -x mm9\_index --threads 8 -U SRR3950566\_1.fastq.gz --reorder | samtools view -Shb - > SRR3950566\_1.bam**
[user@cn4224 ~]$ **bowtie2 -x mm9\_index --threads 8 -U SRR3950566\_2.fastq.gz --reorder | samtools view -Shb - > SRR3950566\_2.bam**
[user@cn4224 ~]$ **bowtie2 -x mm9\_index --threads 8 -U SRR3950559\_1.fastq.gz --reorder | samtools view -Shb - > SRR3950559\_1.bam**
[user@cn4224 ~]$ **bowtie2 -x mm9\_index --threads 8 -U SRR3950559\_2.fastq.gz --reorder | samtools view -Shb - > SRR3950559\_2.bam**
[user@cn4224 ~]$ **bowtie2 -x mm9\_index --threads 8 -U SRR3950560\_1.fastq.gz --reorder | samtools view -Shb - > SRR3950560\_1.bam**
[user@cn4224 ~]$ **bowtie2 -x mm9\_index --threads 8 -U SRR3950560\_2.fastq.gz --reorder | samtools view -Shb - > SRR3950560\_2.bam**

[user@cn4224 ~]$ **hicBuildMatrix --samFiles SRR3950565\_1.bam SRR3950565\_2.bam \
 --binSize 1000 \
 --restrictionSequence GATC \
 --outFileName SRR3950565.cool \
 --QCfolder SRR3950565\_QC \
 --threads 6**
[user@cn4224 ~]$ **hicBuildMatrix --samFiles SRR3950566\_1.bam SRR3950566\_2.bam \
 --binSize 1000 --restrictionSequence GATC \
 --outFileName SRR3950566.cool \
 --QCfolder SRR3950566\_QC \
 --threads 6**
[user@cn4224 ~]$ **hicBuildMatrix --samFiles SRR3950559\_1.bam SRR3950559\_2.bam \
 --binSize 1000 --restrictionSequence GATC \
 --outFileName SRR3950559.cool \
 --QCfolder SRR3950559\_QC \
 --threads 6**
[user@cn4224 ~]$ **hicBuildMatrix --samFiles SRR3950560\_1.bam SRR3950560\_2.bam \
 --binSize 1000 --restrictionSequence GATC \
 --outFileName SRR3950560.cool \
 --QCfolder SRR3950560\_QC \
 --threads 6**
[...]

[user@cn4224 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226

[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. hicexplorer.sh) similar to the following.



```

#! /bin/bash

set -e

module load hicexplorer

cd /lscratch/${SLURM_JOB_ID}

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR395/005/SRR3950565/SRR3950565_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR395/005/SRR3950565/SRR3950565_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR395/006/SRR3950566/SRR3950566_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR395/006/SRR3950566/SRR3950566_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR395/009/SRR3950559/SRR3950559_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR395/009/SRR3950559/SRR3950559_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR395/000/SRR3950560/SRR3950560_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR395/000/SRR3950560/SRR3950560_2.fastq.gz

bowtie2 -x mm9_index --threads 8 -U SRR3950565_1.fastq.gz --reorder | samtools view -Shb - > SRR3950565_1.bam
bowtie2 -x mm9_index --threads 8 -U SRR3950565_2.fastq.gz --reorder | samtools view -Shb - > SRR3950565_2.bam
bowtie2 -x mm9_index --threads 8 -U SRR3950566_1.fastq.gz --reorder | samtools view -Shb - > SRR3950566_1.bam
bowtie2 -x mm9_index --threads 8 -U SRR3950566_2.fastq.gz --reorder | samtools view -Shb - > SRR3950566_2.bam
bowtie2 -x mm9_index --threads 8 -U SRR3950559_1.fastq.gz --reorder | samtools view -Shb - > SRR3950559_1.bam
bowtie2 -x mm9_index --threads 8 -U SRR3950559_2.fastq.gz --reorder | samtools view -Shb - > SRR3950559_2.bam
bowtie2 -x mm9_index --threads 8 -U SRR3950560_1.fastq.gz --reorder | samtools view -Shb - > SRR3950560_1.bam
bowtie2 -x mm9_index --threads 8 -U SRR3950560_2.fastq.gz --reorder | samtools view -Shb - > SRR3950560_2.bam

hicBuildMatrix --samFiles SRR3950565_1.bam SRR3950565_2.bam  \
                                   --binSize 1000 \
                                   --restrictionSequence GATC \
                                   --outFileName SRR3950565.cool \
                                   --QCfolder SRR3950565_QC \
                                   --threads 6
hicBuildMatrix --samFiles SRR3950566_1.bam SRR3950566_2.bam  \
                                   --binSize 1000 --restrictionSequence GATC \
                                   --outFileName SRR3950566.cool \
                                   --QCfolder SRR3950566_QC \
                                   --threads 6
hicBuildMatrix --samFiles SRR3950559_1.bam SRR3950559_2.bam \
                                   --binSize 1000 --restrictionSequence GATC \
                                   --outFileName SRR3950559.cool \
                                   --QCfolder SRR3950559_QC \
                                   --threads 6
hicBuildMatrix --samFiles SRR3950560_1.bam SRR3950560_2.bam \
                                   --binSize 1000 --restrictionSequence GATC \
                                   --outFileName SRR3950560.cool \
                                   --QCfolder SRR3950560_QC \
                                   --threads 6

```

Submit these jobs using the Slurm [sbatch](/docs/userguide.html) command.


Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile for the first step of the pipeline (e.g. hicexplorer.swarm). For example:



```

hicBuildMatrix --samFiles SRR3950565_1.bam SRR3950565_2.bam  \
                                   --binSize 1000 \
                                   --restrictionSequence GATC \
                                   --outFileName SRR3950565.cool \
                                   --QCfolder SRR3950565_QC \
                                   --threads 6
hicBuildMatrix --samFiles SRR3950566_1.bam SRR3950566_2.bam  \
                                   --binSize 1000 --restrictionSequence GATC \
                                   --outFileName SRR3950566.cool \
                                   --QCfolder SRR3950566_QC \
                                   --threads 6
hicBuildMatrix --samFiles SRR3950559_1.bam SRR3950559_2.bam \
                                   --binSize 1000 --restrictionSequence GATC \
                                   --outFileName SRR3950559.cool \
                                   --QCfolder SRR3950559_QC \
                                   --threads 6
hicBuildMatrix --samFiles SRR3950560_1.bam SRR3950560_2.bam \
                                   --binSize 1000 --restrictionSequence GATC \
                                   --outFileName SRR3950560.cool \
                                   --QCfolder SRR3950560_QC \
                                   --threads 6

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f hicexplorer.swarm [-g #] --module hicexplorer
```

where


|  |  |  |  |
| --- | --- | --- | --- |
| -g #  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | --module hicexplorer  Loads the hicexplorer module for each subjob in the swarm
 | |
 | |








