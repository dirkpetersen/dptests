

document.querySelector('title').textContent = 'Canvas on Biowulf';
Canvas on Biowulf


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



Canvas is a tool for calling copy number variants (CNVs) from human DNA sequencing data. It can work either with germline data, or paired tumor/normal samples. Its primary input is aligned reads (in .bam format), and its primary output is a report (in a .vcf file) giving the copy number status of the genome.


Canvas is developed by Illumina. 


Canvas is developed on a Windows platform. It is built as a Singularity container on Biowulf. The script Canvas.sh will start up the Singularity container, and run your command within it, as shown in the following examples. If you follow the examples described on the Illumina site, note that your command on Biowulf should start with the Canvas.dll command -- you do not need 'dotnet' or '/CanvasDIR/'.


### References:


* [Roller, Eric, et al. "Canvas: versatile and scalable detection of copy number variants." *Bioinformatics* 32.15 (2016): 2375-2377.](https://academic.oup.com/bioinformatics/article/32/15/2375/1743834?login=true)


Documentation
* [Canvas Main Site](https://github.com/Illumina/canvas)


Important Notes
* Module Name: Canvas (see [the modules page](/apps/modules.html) for more information)
* Multithreaded application 
* Reference data in /fdb/Canvas/



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive -c 16 --mem 40g**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **mkdir -p /data/$USER/CANVAS\_TEST**

[user@cn3144 ~]$ **cd /data/$USER/CANVAS\_TEST**

[user@cn3144 ~]$ **cp /fdb/platinum\_genomes/strelka\_vcf/variants.vcf.gz .**

[user@cn3144 ~]$ **module load Canvas**

[user@cn3144 ~]$ **Canvas.dll -h**
Usage: Canvas.exe [MODE] [OPTIONS]+
 
Available modes:
	Germline-WGS - CNV calling of a germline sample from whole genome sequencing data
	Somatic-Enrichment - CNV calling of a somatic sample from targeted sequencing data
	Somatic-WGS - CNV calling of a somatic sample from whole genome sequencing data
	Tumor-normal-enrichment - CNV calling of a tumor/normal pair from targeted sequencing data
	SmallPedigree-WGS - CNV calling of a small pedigree from whole genome sequencing data
 
Options:
  -h, --help                 show this message and exit
  -v, --version              print version and exit

[user@cn3144 ~]$ **Canvas.dll SmallPedigree-WGS -b /fdb/platinum\_genomes/bam/sorted.bam \
-r /fdb/Canvas/hg19/Sequence/kmer.fa \
-g /fdb/Canvas/hg19/Sequence/WholeGenomeFasta \
--sample-b-allele-vcf=/data/$USER/CANVAS\_TEST/variants.vcf.gz \
-f /fdb/Canvas/hg19/Sequence/filter13.bed \
-o /data/$USER/CANVAS\_TEST/out**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. submit.sh). For example:



```

#!/bin/bash
set -e
module load Canvas
Canvas.dll Germline-WGS -b /fdb/platinum_genomes/bam/sorted.bam \ 
-r /fdb/Canvas/hg19/Sequence/kmer.fa \
-g /fdb/Canvas/hg19/Sequence/WholeGenomeFasta \
--sample-b-allele-vcf=/data/teacher/canvas/variants.vcf.gz \
-f /fdb/Canvas/hg19/Sequence/filter13.bed \
-o /data/teacher/canvas/out 

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=16 --mem=30g --time 4:00:00 submit.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. canvas.swarm). For example:



```

Canvas.dll Germline-WGS -b sample1.bam --sample-b-allele-vcf=sample1.vcf -o sample1 ...
Canvas.dll Germline-WGS -b sample2.bam --sample-b-allele-vcf=sample2.vcf -o sample2 ...
Canvas.dll Germline-WGS -b sample3.bam --sample-b-allele-vcf sample3.vcf -o sample3 ...

[...] rest of the required options (see batch script above) 

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f canvas.swarm [-g 30] [-t 16] --module Canvas
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module Canvas Loads the Canvas module for each subjob in the swarm 
 | |
 | |
 | |
















