

document.querySelector('title').textContent = 'humann on Biowulf';
humann on Biowulf


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


*From the humann home page:*



>  HUMAnN is a pipeline for efficiently and accurately profiling the
> presence/absence and abundance of microbial pathways in a community from
> metagenomic or metatranscriptomic sequencing data (typically millions of short
> DNA/RNA reads). This process, referred to as functional profiling, aims to
> describe the metabolic potential of a microbial community and its members. More
> generally, functional profiling answers the question "What are the microbes in
> my community-of-interest doing (or capable of doing)?" 



Documentation
* [GitHub](https://github.com/biobakery/humann)
* [Home page](https://huttenhower.sph.harvard.edu/humann)


Important Notes
* Module Name: humann (see [the modules page](/apps/modules.html) 
 for more information)
* For historical reasons there are also modules called humann2 for the series 2 versions
 of humann. They point to the same installation as the equivalent humann modules.
* humann is a multithreaded application
* Reference data in /fdb/humann2/
* Example data in `$HUMANN_TEST_DATA`



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive -c12 --mem=24g --gres=lscratch:100**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load humann/3.6.0**
[user@cn3144 ~]$ **cd /lscratch/$SLURM\_JOB\_ID**
[user@cn3144 ~]$ **cp -r ${HUMANN\_TEST\_DATA:-none} demo**
[user@cn3144 ~]$ **humann --threads $SLURM\_CPUS\_PER\_TASK --input demo/demo.fastq --output demo.out**
[user@cn3144 ~]$ # note that for the humann series 2 the command is `humann2`
WARNING: While bind mounting '/gs10:/gs10': destination is already in the mount point list
Creating output directory: /lscratch/46116226/demo.out
Output files will be written to: /lscratch/46116226/demo.out

Running metaphlan ........

Found g__Bacteroides.s__Bacteroides_dorei : 57.96% of mapped reads
Found g__Bacteroides.s__Bacteroides_vulgatus : 42.04% of mapped reads

Total species selected from prescreen: 2

Selected species explain 100.00% of predicted community composition


Creating custom ChocoPhlAn database ........


Running bowtie2-build ........


Running bowtie2 ........

Total bugs from nucleotide alignment: 2
g__Bacteroides.s__Bacteroides_vulgatus: 1274 hits
g__Bacteroides.s__Bacteroides_dorei: 1318 hits

Total gene families from nucleotide alignment: 548

Unaligned reads after nucleotide alignment: 87.6571428571 %


Running diamond ........


Aligning to reference database: uniref90_201901b_full.dmnd

Total bugs after translated alignment: 3
g__Bacteroides.s__Bacteroides_vulgatus: 1274 hits
g__Bacteroides.s__Bacteroides_dorei: 1318 hits
unclassified: 1599 hits

Total gene families after translated alignment: 815

Unaligned reads after translated alignment: 80.6190476190 %


Computing gene families ...

Computing pathways abundance and coverage ...

Output files created:
/lscratch/46116226/demo.out/demo_genefamilies.tsv
/lscratch/46116226/demo.out/demo_pathabundance.tsv
/lscratch/46116226/demo.out/demo_pathcoverage.tsv

[user@cn3144 ~]$ **ls -lh demo.out**
total 168K
-rw-r--r-- 1 user group 104K Dec 11 12:25 demo_genefamilies.tsv
drwxr-xr-x 2 user group 4.0K Dec 11 12:25 demo_humann_temp
-rw-r--r-- 1 user group 1.4K Dec 11 12:25 demo_pathabundance.tsv
-rw-r--r-- 1 user group 1.3K Dec 11 12:25 demo_pathcoverage.tsv
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. humann.sh), which uses the input file 'humann2.in'. For example:



```

#! /bin/bash

module load humann/3.0.0 || exit 1
cd /lscratch/$SLURM_JOB_ID || exit 1
cp $HUMANN_TEST_DATA/demo.fastq .
mkdir out

# for humann version 2 modules, this command would be humann2
humann --threads $SLURM_CPUS_PER_TASK \
  --input demo.fastq \
  --output out

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] humann.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. humann.swarm). For example:



```

humann --input bample1.bam --output sample1.out
humann --input bample2.bam --output sample2.out
humann --input bample3.bam --output sample3.out

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f humann.swarm -g 10 -t 4 --module humann/3.0.0
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g #  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t #  Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module humann/XXXX  Loads the humann module for each subjob in the swarm 
 | |
 | |
 | |








