

document.querySelector('title').textContent = 'SpeedSeq on Biowulf';
SpeedSeq on Biowulf


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



SpeedSeq is a genome analysis platform designed for rapid whole-genome variant detection and interpretation



### References:


* Chiang C, Layer RM, Faust GG, Lindberg MR, Rose DB, Garrison EP, Marth GT, Quinlan AR, Hall IM. SpeedSeq: ultra-fast personal genome analysis and interpretation. Nature methods. 2015 Aug 10;12(10):966. doi: [10.1038/nmeth.3505](https://doi.org/10.1038/nmeth.3505)


Documentation
* [SpeedSeq Main Site](https://github.com/hall-lab/speedseq)


Important Notes
* Module Name: speedseq (see [the modules page](/apps/modules.html) for more information)
* Environment variables set 
	+ SPEEDSEQ\_HOME* Example files in $SPEEDSEQ\_HOME/example



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**) based on  
 $SPEEDSEQ\_HOME/example/run\_speedseq.sh:



```

[user@biowulf]$ **sinteractive --cpus-per-task 2 --gres lscratch:10**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load speedseq**
[+] Loading VEP 95 on cn3092 
[+] Loading singularity  on cn3092 
[+] Loading ROOT 6.13.02  ... 
[+] Loading gcc  7.3.0  ... 
[+] Loading cnvnator  0.3.3 
[+] Loading speedseq, version 0.1.2-20180208-4e60002... 
[user@cn3144 ~]$ **speedseq align \
 -o example \
 -M $(expr $SLURM\_MEM\_PER\_NODE / 1000) \
 -t $SLURM\_CPUS\_PER\_TASK \
 -T /lscratch/$SLURM\_JOB\_ID/speedseq \
 -p \
 -R "@RG\tID:NA12878\tSM:NA12878\tLB:lib1" \
 $SPEEDSEQ\_HOME/example/data/human\_g1k\_v37\_20\_42220611-42542245.fasta \
 $SPEEDSEQ\_HOME/example/data/NA12878.20slice.30X.fastq.gz**
[user@cn3144 ~]$ **speedseq var \
 -T /lscratch/$SLURM\_JOB\_ID/speedseq \
 -t $SLURM\_CPUS\_PER\_TASK \
 -o example \
 $SPEEDSEQ\_HOME/example/data/human\_g1k\_v37\_20\_42220611-42542245.fasta \
 example.bam**
[user@cn3144 ~]$ **speedseq sv \
 -T /lscratch/$SLURM\_JOB\_ID/speedseq \
 -o example \
 -B example.bam \
 -S example.splitters.bam \
 -D example.discordants.bam \
 -R $SPEEDSEQ\_HOME/example/data/human\_g1k\_v37\_20\_42220611-42542245.fasta**
[user@cn3144 ~]$ **speedseq realign \
 -t $SLURM\_CPUS\_PER\_TASK \
 -T /lscratch/$SLURM\_JOB\_ID/speedseq \
 -o example.realign \
 -M $(expr $SLURM\_MEM\_PER\_NODE / 1000) \
 $SPEEDSEQ\_HOME/example/data/human\_g1k\_v37\_20\_42220611-42542245.fasta \
 example.bam**
[user@cn3144 ~]$ 
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. speedseq.sh). For example:



```

#!/bin/sh
set -e
module load speedseq

test -n "$SLURM_CPUS_PER_TASK" || SLURM_CPUS_PER_TASK=2

speedseq align \
    -o example \
    -M $(expr $SLURM_MEM_PER_NODE / 1000) \
    -t $SLURM_CPUS_PER_TASK \
    -T /lscratch/$SLURM_JOB_ID/speedseq \
    -p \
    -R "@RG\tID:NA12878\tSM:NA12878\tLB:lib1" \
    $SPEEDSEQ_HOME/example/data/human_g1k_v37_20_42220611-42542245.fasta \
    $SPEEDSEQ_HOME/example/data/NA12878.20slice.30X.fastq.gz

speedseq var \
    -T /lscratch/$SLURM_JOB_ID/speedseq \
    -t $SLURM_CPUS_PER_TASK \
    -o example \
    $SPEEDSEQ_HOME/example/data/human_g1k_v37_20_42220611-42542245.fasta \
    example.bam

speedseq sv \
    -T /lscratch/$SLURM_JOB_ID/speedseq \
    -o example \
    -B example.bam \
    -S example.splitters.bam \
    -D example.discordants.bam \
    -R $SPEEDSEQ_HOME/example/data/human_g1k_v37_20_42220611-42542245.fasta

speedseq realign \
    -t $SLURM_CPUS_PER_TASK \
    -T /lscratch/$SLURM_JOB_ID/speedseq \
    -o example.realign \
    -M $(expr $SLURM_MEM_PER_NODE / 1000) \
    $SPEEDSEQ_HOME/example/data/human_g1k_v37_20_42220611-42542245.fasta \
    example.bam

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] --gres lscratch:# speedseq.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. speedseq.swarm). For example:



```

speedseq align \
    -o example1 \
    -M $(expr $SLURM_MEM_PER_NODE / 1000) \
    -t $SLURM_CPUS_PER_TASK \
    -T /lscratch/$SLURM_JOB_ID/speedseq \
    -p \
    -R "@RG\tID:NA12878\tSM:NA12878\tLB:lib1" \
    $SPEEDSEQ_HOME/example/data/human_g1k_v37_20_42220611-42542245.fasta \
    $SPEEDSEQ_HOME/example/data/NA12878.20slice.30X.fastq.gz
speedseq align \
    -o example2 \
    -M $(expr $SLURM_MEM_PER_NODE / 1000) \
    -t $SLURM_CPUS_PER_TASK \
    -T /lscratch/$SLURM_JOB_ID/speedseq \
    -p \
    -R "@RG\tID:NA12878\tSM:NA12878\tLB:lib1" \
    $SPEEDSEQ_HOME/example/data/human_g1k_v37_20_42220611-42542245.fasta \
    $SPEEDSEQ_HOME/example/data/NA12878.20slice.30X.fastq.gz
speedseq align \
    -o example3 \
    -M $(expr $SLURM_MEM_PER_NODE / 1000) \
    -t $SLURM_CPUS_PER_TASK \
    -T /lscratch/$SLURM_JOB_ID/speedseq \
    -p \
    -R "@RG\tID:NA12878\tSM:NA12878\tLB:lib1" \
    $SPEEDSEQ_HOME/example/data/human_g1k_v37_20_42220611-42542245.fasta \
    $SPEEDSEQ_HOME/example/data/NA12878.20slice.30X.fastq.gz
speedseq align \
    -o example4 \
    -M $(expr $SLURM_MEM_PER_NODE / 1000) \
    -t $SLURM_CPUS_PER_TASK \
    -T /lscratch/$SLURM_JOB_ID/speedseq \
    -p \
    -R "@RG\tID:NA12878\tSM:NA12878\tLB:lib1" \
    $SPEEDSEQ_HOME/example/data/human_g1k_v37_20_42220611-42542245.fasta \
    $SPEEDSEQ_HOME/example/data/NA12878.20slice.30X.fastq.gz

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f speedseq.swarm [-g #] -t # --gres lscratch:# --module speedseq
```

where


|  |  |  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --gres lscratch:# | Number of Gigabytes of local scratch space to allocate |
| --module speedseq Loads the TEMPLATE module for each subjob in the swarm 
 | |
 | |
 | |








