

document.querySelector('title').textContent = 'Gridss on Biowulf';
Gridss on Biowulf


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



Gridss is a collection of tools for the detection of genomic rearrangements. It includes a genome-wide break-end assembler, as well as a structural variation caller for Illumina sequencing data. Gridss calls variants based on alignment-guided positional de Bruijn graph genome-wide break-end assembly, split read, and read pair evidence.



References
* Cameron DL, Schröder J, Penington JS, et al. [GRIDSS: sensitive and specific genomic rearrangement detection using positional de Bruijn graph assembly.](https://pubmed.ncbi.nlm.nih.gov/29097403/) , Genome Res. 2017;27(12):2050-2060.


Documentation
* [Gridss github page](https://github.com/PapenfussLab/gridss)


Important Notes
* Module Name: gridss (see [the modules page](/apps/modules.html) for more information). 

* Gridss is installed as a singularity container and contains all the dependencies (e.g., picard, bwa) required for gridss to run.



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) with 5 GB of local disk and run the program. Sample session below::



```

[user@biowulf]$ **sinteractive --mem=8g --cpus-per-task=8 --gres=lscratch:8**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load gridss**
[+] Loading gridss  2.9.4  on cn3144
[+] Loading singularity  3.6.1  on cn3144 

[user@cn3144 ~]$ **cd /lscratch/${SLURM\_JOB\_ID}**

[user@cn3144 ~]$ **cp /fdb/igenomes/Homo\_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa .**

[user@cn3144 ~]$ **cp /data/$USER/test\_dna.bam .**

[user@cn3144 ~]$ **gridss.sh -r genome.fa \
 -o hg38\_test\_output.vcf.gz \
 -a /lscratch/${SLURM\_JOB\_ID}/assembly.bam \
 -w /lscratch/${SLURM\_JOB\_ID} \
 test\_dna.bam** 

Using working directory "/lscratch/46116226"
Wed Aug 19 13:40:00 EDT 2020: Full log file is: gridss.full.20200819_134000.cn3144.5147.log
Wed Aug 19 13:40:00 EDT 2020: Found /usr/bin/time
Wed Aug 19 13:40:00 EDT 2020: Using reference genome "genome.fa"
Wed Aug 19 13:40:00 EDT 2020: Using assembly bam /lscratch/46116226
Wed Aug 19 13:40:00 EDT 2020: Using output VCF hg38_test_output
Wed Aug 19 13:40:00 EDT 2020: Using 8 worker threads.
Wed Aug 19 13:40:00 EDT 2020: Using no blacklist bed. The encode DAC blacklist is recommended for hg19.
Wed Aug 19 13:40:00 EDT 2020: Using JVM maximum heap size of 25g for assembly and variant calling.
Wed Aug 19 13:40:00 EDT 2020: Using input file test_dna.bam
Wed Aug 19 13:40:00 EDT 2020: Found /usr/bin/Rscript
Wed Aug 19 13:40:00 EDT 2020: Found /usr/bin/samtools
Wed Aug 19 13:40:00 EDT 2020: Found /usr/bin/java
Wed Aug 19 13:40:00 EDT 2020: Found /usr/bin/bwa
Wed Aug 19 13:40:01 EDT 2020: samtools version: 1.7+htslib-1.7-2
Wed Aug 19 13:40:01 EDT 2020: R version: R scripting front-end version 3.6.3 (2020-02-29)
Wed Aug 19 13:40:01 EDT 2020: bwa Version: 0.7.17-r1188

[...]

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. gridss.sh). For example:



```

#!/bin/bash
set -e
module load gridss

gridss.sh -r genome.fa \
          -o hg38_test_output.vcf.gz \
          -a /lscratch/${SLURM_JOB_ID}/assembly.bam \
          -w /lscratch/${SLURM_JOB_ID} \
          test_dna.bam 


```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch  [--mem=#] --gres=lscratch:5 gridss.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.


Create a swarmfile (e.g. gridss.swarm). For example:



```

gridss.sh -r genome.fa -o hg38_test_output1.vcf.gz -a /lscratch/${SLURM_JOB_ID}/assembly1.bam \
          -w /lscratch/${SLURM_JOB_ID} test_dna1.bam
gridss.sh -r genome.fa -o hg38_test_output2.vcf.gz -a /lscratch/${SLURM_JOB_ID}/assembly2.bam \
           -w /lscratch/${SLURM_JOB_ID} test_dna2.bam
gridss.sh -r genome.fa -o hg38_test_output3.vcf.gz -a /lscratch/${SLURM_JOB_ID}/assembly3.bam \
           -w /lscratch/${SLURM_JOB_ID} test_dna3.bam

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f gridss.swarm [-t #] --gres=lscratch:5 --module gridss
```

where


|  |  |  |  |
| --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | --gres=lscratch:5 allocate 5 GB of local disk for each swarm subjob
 | |
 | |








