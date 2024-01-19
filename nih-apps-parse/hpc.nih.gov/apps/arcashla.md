

document.querySelector('title').textContent = 'arcasHLA on HPC';
arcasHLA on HPC


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

  arcasHLA: high resolution HLA typing from RNA seq. arcasHLA performs high resolution genotyping for HLA class I and class II genes from RNA sequencing, supporting both paired and single-end samples.

 ### 


Documentation * <https://github.com/RabadanLab/arcasHLA>



Important Notes * Module Name: arcashla (see [the modules 
 page](/apps/modules.html) for more information)
* example : /usr/local/apps/arcashla/0.5.0/test





Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=5g --cpus-per-task=4**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load arcashla**
[user@cn3144 ~]$ **cd /data/$USER/**
[user@cn3144 ~]$ **cp -r /usr/local/apps/arcashla/0.5.0/test .**
[user@cn3144 ~]$ **arcasHLA extract test/test.bam -o test/output -t $SLURM\_CPUS\_PER\_TASK**
[user@cn3144 ~]$ **arcasHLA genotype test/output/test.extracted.1.fq.gz test/output/test.extracted.2.fq.gz -g A,B,C,DPB1,DQB1,DQA1,DRB1 -o test/output -t $SLURM\_CPUS\_PER\_TASK**
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```




Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. batch.sh). For example:



```

#!/bin/bash
set -e
module load arcashla
cd /data/$USER/test
arcasHLA extract test.bam -o output -t $SLURM_CPUS_PER_TASK
arcasHLA genotype output/test.extracted.1.fq.gz output/test.extracted.2.fq.gz -g A,B,C,DPB1,DQB1,DQA1,DRB1 -o output -t $SLURM_CPUS_PER_TASK
```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--mem=#] [--cpus-per-task=#] batch.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. job.swarm). For example:



```

cd dir1; arcasHLA extract test.bam -o output -t $SLURM_CPUS_PER_TASK
cd dir2; arcasHLA extract test.bam -o output -t $SLURM_CPUS_PER_TASK
cd dir3; arcasHLA extract test.bam -o output -t $SLURM_CPUS_PER_TASK

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f job.swarm [-g #] [-t #] --module arcashla
```

where
 

|  |  |
| --- | --- |
| -g *#*  | Number of Gigabytes of memory required for each process 
 (1 line in the swarm command file)  |
| -t *#*  | Number of threads required for each job
 (1 line in the swarm command file)  |
| --module  | Loads the module for each subjob in the swarm  |




