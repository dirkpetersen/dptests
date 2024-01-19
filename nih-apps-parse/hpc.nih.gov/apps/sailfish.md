

document.querySelector('title').textContent = 'sailfish on Biowulf';
sailfish on Biowulf


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



Sailfish quantifies the expression of a given set of transcripts using NGS reads. 
It is run in two stages: (1) The indexing step is run once per set of transcripts
(2) The quantification step is run once for each sample.



### References:


* R. Patro, S. M. Mount, and C. Kingsford. *Sailfish enables 
 alignment-free isoform quantification from RNA-seq reads using 
 lightweight algorithms*. Nature Biotechnology 2014, 32:462-464.
 [Pubmed](http://www.ncbi.nlm.nih.gov/pubmed/24752080) | 
 [PMC](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4077321/) | 
 [Journal](http://www.nature.com/nbt/journal/v32/n5/full/nbt.2862.html)


Documentation
* [GitHub](https://github.com/kingsfordgroup/sailfish)
* [Manual](http://sailfish.readthedocs.io/en/master/index.html)


Important Notes
* Module Name: sailfish (see [the modules page](/apps/modules.html) for more information)
* sailfish is a multithreaded application. Make sure to match the number of cpus requested with
 the number of threads (-p).
* Example files in `SAILFISH_TEST_DATA`



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=8g --cpus-per-task=4 --gres=lscratch:20**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144]$ **cd /lscratch/$SLURM\_JOB\_ID**
[user@cn3144]$ **module load sailfish**
[user@cn3144]$ **zcat $SAILFISH\_TEST\_DATA/gencode.vM9.transcripts.fa.gz > M9.fa**
[user@cn3144]$ **sailfish index -t M9.fa -o M9.idx -p $SLURM\_CPUS\_PER\_TASK**
[user@cn3144]$ **cp $SAILFISH\_TEST\_DATA/ENCFF138LJO.fastq.gz .**
[user@cn3144]$ **sailfish quant -i M9.idx -r <(zcat ENCFF138LJO.fastq.gz) --libType U \
 -o quant -p $SLURM\_CPUS\_PER\_TASK**

[user@cn3144]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. sailfish.sh), which uses the input file 'sailfish.in'. For example:



```

#! /bin/bash

module load sailfish/0.10.0 || exit 1
wd=$PWD
cd /lscratch/$SLURM_JOB_ID || exit 1

# get transcriptome from the example directory
# This is usually done once - not for each job. Only included
# here to show all steps involved in sailfish quantitation.
zcat $SAILFISH_TEST_DATA/gencode.vM9.transcripts.fa.gz \
    > gencode.vM9.transcripts.fa

# index the transcripts
sailfish index -t gencode.vM9.transcripts.fa -o gencode.vM9.idx \
    -p $SLURM_CPUS_PER_TASK

# quantify the transcripts
cp $SAILFISH_TEST_DATA/ENCFF138LJO.fastq.gz .
sailfish quant -i gencode.vM9.idx -l U \
    -r <(zcat ENCFF138LJO.fastq.gz) \
    -o quant -p $SLURM_CPUS_PER_TASK
cp -r quant $wd

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=8 --mem=8g --gres=lscratch:16 sailfish.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. sailfish.swarm). For example:



```

sailfish quant -i gencode.vM9.idx -l U -r <(zcat sample1.fq.gz) \
    -o quant_sample1 -p $SLURM_CPUS_PER_TASK
sailfish quant -i gencode.vM9.idx -l U -r <(zcat sample2.fq.gz) \
    -o quant_sample2 -p $SLURM_CPUS_PER_TASK
sailfish quant -i gencode.vM9.idx -l U -r <(zcat sample3.fq.gz) \
    -o quant_sample3 -p $SLURM_CPUS_PER_TASK

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f sailfish.swarm -g 8 -t 8 --module sailfish/0.10.0
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g #  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t #  Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module sailfish  Loads the sailfish module for each subjob in the swarm 
 | |
 | |
 | |








