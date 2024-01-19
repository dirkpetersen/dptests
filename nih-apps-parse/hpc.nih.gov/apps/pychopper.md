

document.querySelector('title').textContent = 'pychopper on Biowulf';
pychopper on Biowulf


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


Pychopper is used to identify, orient and trim full-length Nanopore
cDNA reads. The tool is also able to rescue fused reads.



Documentation
* pychopper on [GitHub](https://github.com/nanoporetech/pychopper)


Important Notes
* Module Name: pychopper (see [the modules page](/apps/modules.html) for more information)
* pychopper can use multiple CPUs. Please match your allocation with the number
 of threads used. pychopper does not scale efficiently to more than 8 CPUs
* pychopper ≤ 2.4.0 cannot natively read gzip input. ≥2.7.1 can read gzip files but
 it is inefficient. It remains more efficient to uncompress input to lscratch for processing.
* pychopper often reads the fastq file twice. Therefore it's best to move/unpack fastq files
 into lscratch. See example below.
* `-m edlib` is associated with stalled runs and should probably be avoided.
* Example files in `$PYCHOPPER_TEST_DATA`
* the command `cdna_classifier.py` was renamed to `pychopper` between 2.4.0 and 2.7.1



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=10g --cpus-per-task=6 --gres=lscratch:10**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144]$ **cd /lscratch/$SLURM\_JOB\_ID**
[user@cn3144]$ **module load pychopper**
[user@cn3144]$ **zcat $PYCHOPPER\_TEST\_DATA/SIRV\_E0\_pcs109\_25k.fq.gz > input.fastq**
[user@cn3144]$ ## pychopper was called cdna_classifier.py in versions ≤ 2.4.0
[user@cn3144]$ **pychopper -r report.pdf -u unclassified.fastq -t $SLURM\_CPUS\_PER\_TASK \
 -w rescued.fastq input.fastq - | gzip -c > /data/$USER/temp/full\_length.fastq.gz**
Using kit: PCS109
Configurations to consider: "+:SSP,-VNP|-:VNP,-SSP"
Total fastq records in input file: 25000
Tuning the cutoff parameter (q) on 9465 sampled reads (40.0%) passing quality filters (Q ≥ 7.0).
Optimizing over 30 cutoff values.
100%|████████████████████████████████████████████████████████| 30/30
Best cutoff (q) value is 0.3448 with 88% of the reads classified.
Processing the whole dataset using a batch size of 4166:
 94%|██████████████████████████████████████████████████      | 23614/25000
Finished processing file: input.fastq
Input reads failing mean quality filter (Q < 7.0): 1386 (5.54%)
Output fragments failing length filter (length < 50): 0
-----------------------------------
Reads with two primers: 86.93%
Rescued reads:          3.16%
Unusable reads:         9.91%
-----------------------------------


```

Move the rescuted and unclassified reads and the reports if you need them before
ending the session.



```

[user@cn3144]$ **gzip -c rescued.fastq > /data/$USER/temp/rescued.fastq.gz**
[user@cn3144]$ **gzip -c unclassified.fastq > /data/$USER/temp/unclassified.fastq.gz**
[user@cn3144]$ **mv report.pdf /data/$USER/temp**
[user@cn3144]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. pychopper.sh), which uses the input file 'pychopper.in'. For example:



```

#!/bin/bash
module load pychopper/2.7.1 || exit 1
cd /lscratch/$SLURM_JOB_ID
zcat $PYCHOPPER_TEST_DATA/SIRV_E0_pcs109_25k.fq.gz > input.fastq
## use cdna_classifier.py instead of pychopper for versions ≤ 2.4.0
pychopper -r report.pdf -u unclassified.fastq -t $SLURM_CPUS_PER_TASK \
    -w rescued.fastq input.fastq - | gzip -c > /data/$USER/temp/full_length.fastq.gz
gzip -c rescued.fastq > /data/$USER/temp/rescued.fastq.gz
gzip -c unclassified.fastq > /data/$USER/temp/unclassified.fastq.gz
mv report.pdf /data/$USER/temp

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=6 --mem=10g pychopper.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. pychopper.swarm). For example:



```

zcat input1.fastq.gz > /lscratch/$SLURM_JOB_ID/input1.fastq && \
  pychopper -t $SLURM_CPUS_PER_TASK /lscratch/$SLURM_JOB_ID/input1.fastq - \
  | gzip -c > /data/$USER/temp/full_length1.fastq.gz
zcat input2.fastq.gz > /lscratch/$SLURM_JOB_ID/input1.fastq && \
  pychopper -t $SLURM_CPUS_PER_TASK /lscratch/$SLURM_JOB_ID/input2.fastq - \
  | gzip -c > /data/$USER/temp/full_length2.fastq.gz

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f pychopper.swarm -g 10 -t 6 --module pychopper/2.0.3
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g #  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t #  Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module pychopper  Loads the pychopper module for each subjob in the swarm 
 | |
 | |
 | |








