

document.querySelector('title').textContent = 'cgpBattenberg on Biowulf';
cgpBattenberg on Biowulf


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


 An installation helper, perl wrapper and the R program Battenberg which
detects subclonality and copy number in matched NGS data.


Documentation
* cgpBattenberg [GitHub repo](https://github.com/cancerit/cgpBattenberg)


Important Notes
* Module Name: cgpBattenberg (see [the modules page](/apps/modules.html) for more information)
* cgpBattenberg.pl can use multiple threads (`-t`). Please match the number of threads with
 * cgpBattenberg is part of the Cancer Genome Project and is closely related to the programs [BRASS](/apps/BRASS.html) and [Ascat NGS](/apps/ascatNgs.html) as well as the utilites [VAGrENT](https://github.com/cancerit/VAGrENT) and [PCAP-core](https://github.com/cancerit/PCAP-core). All of these programs can be added to your path using the cancerit-wgs module. To get the most recent versions of all of these, use the cancerit-wgs/latest module version.
 the number of allocated CPUs
* Reference data in `/fdb/cancerit-wgs/cgpBattenberg/`



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session. Note
that the cgpBattenberg.pl is part of the cancerit-wgs tools.



```

[user@biowulf]$ **sinteractive --cpus-per-task=4**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144]$ **module load cgpBattenberg**
Usage:
    battenberg.pl [options]

      Required parameters:
        -outdir                -o   Folder to output result to.
        -reference             -r   Path to reference genome index file *.fai
        -tumbam                -tb  Path to tumour bam file
                                     - when '-a' defined sample name
        -normbam               -nb  Path to normal bam file
                                     - when '-a' defined sample name
        -gender                -ge  Gender, XX, XY or L (see -gl)
        -impute-info           -e   Location of the impute info file
        -thousand-genomes-loc  -u   Location of the directory containing 1k genomes data
        -ignore-contigs-file   -ig  File containing contigs to ignore
                                    - specifically male sex chromosome, mitochondria and non primary contigs
        -gc-correction-loc     -gc  Path to gc correction files

      Optional parameters:
[...snip...]

[user@cn3144]$ **battenberg.pl -p output \
 -r /fdb/igenomes/Homo\_sapiens/Ensembl/GRCh37/Sequence/WholeGenomeFasta/genome.fa \
 -tb tumor.bam \
 -nb normal.bam \
 -ge XX \
 -impute-info /fdb/cancerit-wgs/cgpBattenberg/impute/impute\_info.txt \
 -thousand-genomes-loc /fdb/cancerit-wgs/cgpBattenberg/1000genomesloci \
 -ignore-contigs-file ignore\_contigs \
 -gc-correction-loc /fdb/cancerit-wgs/cgpBattenberg/battenberg\_wgs\_gc\_correction\_1000g\_v3 \
 -species Human -assembly 37 \
 -t $SLURM\_CPUS\_PER\_TASK**

[user@cn3144]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf]$

```

Notes:


* The ignore contigs file should include the Y, MT, and non-primary contigs. Running without
 this will result in errors.
* Increasing the number of threads results in a proportional increase in the amount of memory
 required. In one test I was able to run 6 threads with 40GB of memory with tumor and normal files
 of approximately 25GB.



Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. cgpBattenberg.sh), which uses the input file 'cgpBattenberg.in'. For example:



```

#!/bin/bash
module load cgpBattenberg
battenberg.pl -p output \
    -r /fdb/igenomes/Homo_sapiens/Ensembl/GRCh37/Sequence/WholeGenomeFasta/genome.fa \
    -tb tumor.bam \
    -nb normal.bam \
    -ge XX \
    -impute-info /fdb/cancerit-wgs/cgpBattenberg/impute/impute_info.txt \
    -thousand-genomes-loc /fdb/cancerit-wgs/cgpBattenberg/1000genomesloci \
    -ignore-contigs-file ignore_contigs \
    -gc-correction-loc /fdb/cancerit-wgs/cgpBattenberg/battenberg_wgs_gc_correction_1000g_v3 \
    -species Human -assembly 37 \
    -t $SLURM_CPUS_PER_TASK

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=16 --mem=30g cgpBattenberg.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. cgpBattenberg.swarm). For example:



```

battenberg.pl -p output -tb tumor1.bam -nb normal1.bam ... -t $SLURM_CPUS_PER_TASK 
battenberg.pl -p output -tb tumor2.bam -nb normal2.bam ... -t $SLURM_CPUS_PER_TASK 
battenberg.pl -p output -tb tumor3.bam -nb normal3.bam ... -t $SLURM_CPUS_PER_TASK 

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f cgpBattenberg.swarm -g 30 -t 16 --module cgpBattenberg
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g #  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t #  Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module cgpBattenberg  Loads the cgpBattenberg module for each subjob in the swarm 
 | |
 | |
 | |








