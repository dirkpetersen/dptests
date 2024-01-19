

document.querySelector('title').textContent = "dorado";
dorado on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
 |



Dorado is a basecaller for Oxford Nanopore reads.




Documentation
* dorado on [GitHub](https://github.com/nanoporetech/dorado)


Important Notes
* Module Name: dorado (see [the modules page](/apps/modules.html) for more information)
* Requires a V100/V100x or newer GPU for basecalling. Alignment is not accelerated.
* Use [pod5](/apps/pod5.html) input format for optimal performance
* Models are found in ${DORADO\_MODELS}
* Example files in ${DORADO\_TEST\_DATA}



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive --gres=gpu:v100:1,lscratch:200 --mem=16g --cpus-per-task=6**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load dorado/0.3.4**
[user@cn3144 ~]$ **cd /lscratch/$SLURM\_JOB\_ID**
[user@cn3144 ~]$ **cp -rL "${DORADO\_TEST\_DATA:-none}" input**
[user@cn3144 ~]$ **ls -lh input**
-rw-r--r--. 1 user group 20G Jun  2 17:21 reads.pod5
[user@cn3144 ~]$ # emits unaligned bam by default
[user@cn3144 ~]$ **dorado basecaller --device cuda:all ${DORADO\_MODELS}/dna\_r9.4.1\_e8\_sup@v3.3 input > output.bam**
[user@cn3144 ~]$ **ls -lh output.bam**
-rw-r--r-- 1 user group 2.1G Jun  2 20:02 output.bam
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. dorado.sh). For example:



```

#!/bin/bash
set -e
module load dorado/0.3.4
cd /lscratch/$SLURM_JOB_ID
cp -rL "${DORADO_TEST_DATA:-none}" input
dorado basecaller --device cuda:all ${DORADO_MODELS}/dna_r9.4.1_e8_sup@v3.3 input > output.bam

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=6 --mem=16g --gres=lscratch:50,gpu:v100x:1 dorado.sh
```







