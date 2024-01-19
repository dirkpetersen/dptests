

document.querySelector('title').textContent = 'Slamdunk on Biowulf';
Slamdunk on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
 |



Streamlining SLAMseq analysis with ultra-high sensitivity: SlamDunk is a novel, fully automated software tool for automated, robust, scalable and reproducible SLAMseq data analysis.



### References:


* [Quantification of experimentally induced nucleotide conversions in high-throughput sequencing datasets](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2849-7), Neumann, T., Herzog, V. A. et al. BMC Bioinformatics, 20(1), 258, 2019.


Documentation
* [Slamdunk documentation](http://t-neumann.github.io/slamdunk/)


Important Notes
* Module Name: slamdunk (see [the modules page](/apps/modules.html) for more information)
* Multithreaded



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. 

  
Sample session (user input in **bold**):




```

[user@biowulf]$ **sinteractive --cpus-per-task=8 --mem=20g --gres=lscratch:20**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load slamdunk**

[user@cn3144 ~]$ **cd /lscratch/$SLURM\_JOB\_ID**

[user@cn3144 ~]$ **cp /fdb/genome/hg19/chr\_all.fa .**

[user@cn3144 ~]$ **slamdunk all -r chr\_all.fa -b ./exampleBAM.bed \
 -o out -t $SLURM\_CPUS\_PER\_TASK ./exampleBAM.bam**
slamdunk all
Running slamDunk map for 1 files (8 threads)
.
Running slamDunk sam2bam for 1 files (8 threads)
.
Creating output directory: out/filter
Running slamDunk filter for 1 files (8 threads)
.
[...]

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. slamdunk.sh). For example:



```

#!/bin/bash
set -e
module load slamdunk
cd /lscratch/$SLURM_JOB_ID
cp /fdb/genome/hg19/chr_all.fa .
slamdunk all -r chr_all.fa -b ./exampleBAM.bed \
             -o out -t $SLURM_CPUS_PER_TASK ./exampleBAM.bam

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] [--gres=lscratch:#] slamdunk.sh
```







