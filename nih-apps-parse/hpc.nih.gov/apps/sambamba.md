

document.querySelector('title').textContent = 'Sambamba on Biowulf';
Sambamba on Biowulf


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


 Sambamba is a high performance modern robust and fast tool (and library), 
 written in the D programming language, for working with SAM and BAM files. 
 Current parallelised functionality is an important subset of samtools functionality, 
 including view, index, sort, markdup, and depth.


### References:


* <https://academic.oup.com/bioinformatics/article/31/12/2032/214758>


Documentation
* <https://github.com/biod/sambamba>



Important Notes
* Module Name: sambamba (see [the 
 modules page](/apps/modules.html) for more information)
* Multithreaded





Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --cpus-per-task=4 --gres=lscratch:10**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load sambamba**
[user@cn3144 ~]$ **mkdir -p /data/$USER/sambamba && cd /data/$USER/sambamba**
[user@cn3144 ~]$ **cp $SAMBAMBA\_TEST\_DATA/\* .**
[user@cn3144 ~]$ **sambamba sort -t $SLURM\_CPUS\_PER\_TASK --tmpdir=lscratch/$SLURM\_JOB\_ID issue\_204.bam**
[user@cn3144 ~]$ **sambamba view -t $SLURM\_CPUS\_PER\_TASK -c -F "proper\_pair" issue\_204.bam**
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```




Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. sambamba.sh). For example:



```

#!/bin/bash
set -e
module load sambamba
sambamba sort -t $SLURM_CPUS_PER_TASK --tmpdir=lscratch/$SLURM_JOB_ID issue_204.bam
```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=4 --gres=lscratch:10 sambamba.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources. Create a swarmfile (e.g. sambamba.swarm). For example:



```
cd dir1; sambamba view --reference-info issue_204.bam
cd dir2; sambamba view --header --format=json issue_204.bam
cd dir3; sambamba flagstat -t $SLURM_CPUS_PER_TASK --show-progress issue_204.bam 
```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f sambamba.swarm -t 4 --module sambamba
```

where
 

|  |  |
| --- | --- |
| -t *#* | Number of threads/CPUs required for each process (1 line 
 in the swarm command file).  |
| --module  | Loads the module for each subjob in the swarm  |




