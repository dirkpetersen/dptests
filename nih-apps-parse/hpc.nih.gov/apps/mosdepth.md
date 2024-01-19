

document.querySelector('title').textContent = 'mosdepth on Biowulf';
mosdepth on Biowulf


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



Fast BAM/CRAM depth calculation for WGS, exome, or targeted sequencing.



### References:


* B. S. Pedersen and A. R. Quinlan. 
 *mosdepth: quick coverage calculation for genomes and exomes.*
 Bioinformatics 2017.
 [PubMed](https://www.ncbi.nlm.nih.gov/pubmed/29096012)  | 
 PMC  | 
 [Journal](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btx699/4583630)


Documentation
* mosdepth Main Site: [GitHub/mosdepth](https://github.com/brentp/mosdepth)


Important Notes
* Module Name: mosdepth (see [the modules page](/apps/modules.html) for more information)
* mosdepth can use multiple threads for bam decompression but does not scale well 
 to more than 4 threads.
* Test data can be found in `$MOSDEPTH_TEST_DATA`



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and follow
the example in the following example session.



```

[user@biowulf]$ **sinteractive --mem=10g --cpus-per-task=4 --gres=lscratch:50**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load mosdepth**
[user@cn3144 ~]$ **cd /lscratch/$SLURM\_JOB\_ID**
[user@cn3144 ~]$ **cp -rL ${MOSDEPTH\_TEST\_DATA:-none} data** # ~20 GB of data
[user@cn3144 ~]$ # whole exome data analyzed by capture probes
[user@cn3144 ~]$ **mosdepth --by data/Agilent\_SureSelect\_Human\_AllExon\_V4\_Covered.bed \
 HG00096\_WES data/HG00096.mapped.ILLUMINA.bwa.GBR.exome.20120522.bam**

[user@cn3144 ~]$ **ls -lh**
-rw-r--r-- 1 user group 336K May  8 18:41 HG00096_WES.mosdepth.global.dist.txt
-rw-r--r-- 1 user group 149K May  8 18:41 HG00096_WES.mosdepth.region.dist.txt
-rw-r--r-- 1 user group 570M May  8 18:41 HG00096_WES.per-base.bed.gz
-rw-r--r-- 1 user group  87K May  8 18:41 HG00096_WES.per-base.bed.gz.csi
-rw-r--r-- 1 user group 2.2M May  8 18:41 HG00096_WES.regions.bed.gz
-rw-r--r-- 1 user group 350K May  8 18:41 HG00096_WES.regions.bed.gz.csi

```

Note that version 0.2.0 produced only a single distribution file (
`HG00096_WES.mosdepth.dist.txt`.



```

[user@cn3144 ~]$ # create some graphs. Output file is dist.html
[user@cn3144 ~]$ **plot-dist.py HG00096\_WES.mosdepth.region.dist.txt**

[user@cn3144 ~]$ # WGS data in 500 nt windows. No per-base data is written (-n)
[user@cn3144 ~]$ **mosdepth -t3 -n --by 500 HG00096\_WGS \
 data/HG00096.mapped.ILLUMINA.bwa.GBR.low\_coverage.20120522.bam** 
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. mosdepth.sh), which uses the input file 'mosdepth.in'. For example:



```

#! /bin/bash

module load mosdepth/0.3.3 || exit 1

cd /lscratch/$SLURM_JOB_ID || exit 1
mkdir data
cp ${MOSDEPTH_TEST_DATA}/Agilent_SureSelect_Human_AllExon_V4_Covered.bed data
cp ${MOSDEPTH_TEST_DATA}/HG00096.mapped.ILLUMINA.bwa.GBR.exome.20120522.bam data

mosdepth -n -t2 --by data/Agilent_SureSelect_Human_AllExon_V4_Covered.bed \ 
    HG00096_WES data/HG00096.mapped.ILLUMINA.bwa.GBR.exome.20120522.bam
plot-dist.py HG00096_WES.mosdepth.region.dist.txt

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --mem=5g mosdepth.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. mosdepth.swarm). For example:



```

mosdepth -t2 -n --by 500 sample1 sample1.bam
mosdepth -t2 -n --by 500 sample2 sample2.bam
mosdepth -t2 -n --by 500 sample3 sample3.bam

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f mosdepth.swarm -g 5 --module mosdepth
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g #  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t #  Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module mosdepth  Loads the mosdepth module for each subjob in the swarm 
 | |
 | |
 | |








