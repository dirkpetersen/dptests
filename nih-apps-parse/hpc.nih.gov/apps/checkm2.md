

document.querySelector('title').textContent = "checkm2";
checkm2 on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
 |



Rapid assessment of genome bin quality using machine learning. From the documentation:




> 
>  Unlike CheckM1, CheckM2 has universally trained machine learning models it
>  applies regardless of taxonomic lineage to predict the completeness and
>  contamination of genomic bins. This allows it to incorporate many lineages in
>  its training set that have few - or even just one - high-quality genomic
>  representatives, by putting it in the context of all other organisms in the
>  training set. As a result of this machine learning framework, CheckM2 is also
>  highly accurate on organisms with reduced genomes or unusual biology, such as
>  the Nanoarchaeota or Patescibacteria.
> 



Documentation
* CheckM2 on [GitHub](https://github.com/chklovski/CheckM2)


Important Notes
* Module Name: checkm2 (see [the modules page](/apps/modules.html) for more information)
* This application is multithreaded. Please match the number of allocated CPUs to the number of threads.
* Example files in $CHECKM2\_TEST\_DATA
* CheckM2 will automatically use lscratch as a temp dir if lscratch has been allocated.
* Benchmarking suggests that CheckM2 may not scale efficiently to more than 16 CPUs



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive --mem=12g --cpus-per-task=16 --gres=lscratch:20**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load checkm2**
[user@cn3144 ~]$ **cd /lscratch/$SLURM\_JOB\_ID**
[user@cn3144 ~]$ **mkdir tmp**
[user@cn3144 ~]$ **cp -Lr $CHECKM2\_TEST\_DATA/fasta .**
[user@cn3144 ~]$ **checkm2 predict --threads=$SLURM\_CPUS\_PER\_TASK \
 --input ./fasta \
 --output-directory ./checkm2**
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```

Timing and efficiency with different numbers of CPUs for the 36 genomic bins in the test
data set:





| CPUs | Memory [GiB] | Runtime [minutes] | Est. efficiency |
| --- | --- | --- | --- |
| 2 | 9 | 29.9 | 100% |
| 4 | 9 | 15.3 | 98% |
| 8 | 9 | 10.5 | 72% |
| 16 | 9 | 6.2 | 60% |




Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. checkm2.sh). For example:



```

#!/bin/bash
set -e
module load checkm2/1.0.2
# uses lscratch automatically as tmp dir
checkm2 predict \
    --threads=$SLURM_CPUS_PER_TASK \
    --input $CHECKM2_TEST_DATA/fasta \
    --output-directory test-checkm2.out

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=8 --mem=12g --gres=lscratch:10 checkm2.sh
```







