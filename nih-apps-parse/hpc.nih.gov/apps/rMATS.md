

document.querySelector('title').textContent = 'rMATS on Biowulf';
rMATS on Biowulf


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



MATS is a computational tool to detect differential alternative splicing events from RNA-Seq data. 
The statistical model of MATS calculates the P-value and false discovery rate 
that the difference in the isoform ratio of a gene between two conditions 
exceeds a given user-defined threshold. 
The replicate MATS (rMATS) is designed for detection of differential alternative splicing 
from replicate RNA-Seq data.



### References:


* Shen S., Park JW., Lu ZX., Lin L., Henry MD., Wu YN., Zhou Q., Xing Y. rMATS: Robust and Flexible Detection of Differential Alternative Splicing from Replicate RNA-Seq Data. PNAS, 111(51):E5593-601. doi: 10.1073/pnas.1419161111


Documentation
* [rMATS home page](http://rnaseq-mats.sourceforge.net/index.html)


Important Notes
* Module Name: rmats (see [the modules page](/apps/modules.html) for more information)
* Multithreaded app (set --nthread $SLURM\_CPUS\_PER\_TASK)
* Example files in /usr/local/apps/rmats/TEST\_DATA* Reference data in /fdb/STAR\_indices* Unusual environment variables set
	+ **RMATS\_HOME**  installation directory
	+ **RMATS\_BIN**  executable directory
	+ **RMATS\_SRC**  source code directory
	+ **RMATS\_DATA**  sample data directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive -c 16 --mem 45g --gres=lscratch:20**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **mkdir -p /data/$USER/rmats && cd /data/$USER/rmats**

```

Here is how one can use the most recent rMATS version 4.1.2:

```

[user@cn3144 ~]$ **module load rmats/4.1.2**
[+] Loading singularity  3.8.5-1  on cn3144
[+] Loading rMATS  4.1.2
[user@cn3144 ~]$  **cp -r $RMATS\_DATA/\* .** 

[user@cn3144 ~]$  **rmats.py --s1 $PWD/s1.txt --s2 $PWD/s2.txt --gtf gtf/Homo\_sapiens.Ensembl.GRCh37.72.gtf --bi /fdb/STAR\_indices/2.7.8a/GENCODE/Gencode\_human/release\_27/genes-100 --od out\_test -t paired --nthread $SLURM\_CPUS\_PER\_TASK --readLength 50 --tophatAnchor 8 --cstat 0.0001 --tstat 6 --tmp /lscratch/${SLURM\_JOB\_ID}** 
...

```

In order to use an older rMATS version 4.0.2:

```

[user@cn3144 ~]$ **module load rmats/4.0.2**
[+] Loading rmats  4.0.2  on cn3144
[+] Loading singularity  3.8.5-1  on cn3144
[user@cn3144 ~]$ **cp -r /usr/local/apps/rmats/testData/\* .**

[user@cn3144 ~]$ **export TMPDIR=/lscratch/$SLURM\_JOBID** # write temp files in lscratch

[user@cn3144 ~]$ **rmats --s1 $PWD/s1.txt --s2 $PWD/s2.txt --gtf gtf/Homo\_sapiens.Ensembl.GRCh37.72.gtf --bi /fdb/STAR\_indices/2.6.1c/GENCODE/Gencode\_human/release\_27/genes-100 --od out\_test -t paired --nthread $SLURM\_CPUS\_PER\_TASK --readLength 50 --tophatAnchor 8 --cstat 0.0001 --tstat 6** 

```


```

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. submit.sh). For example:



```

#!/bin/bash
set -e
module load rmats
export TMPDIR=/lscratch/$SLURM_JOBID
rmats.py  --s1 $PWD/s1.txt --s2 $PWD/s2.txt --gtf $PWD/gtf/Homo_sapiens.Ensembl.GRCh37.72.gtf --bi /fdb/STAR_indices/2.7.8a/GENCODE/Gencode_human/release_27/genes-100 --od out_test -t paired --nthread $SLURM_CPUS_PER_TASK --readLength 50 --tophatAnchor 8 --cstat 0.0001 --tstat 6 --tmp /lscratch/${SLURM_JOB_ID}
```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=16 --mem=30g --gres=lscratch:20 submit.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. rmats.swarm). For example:



```

rmats.py  --s1 $PWD/s1.txt --s2 $PWD/s2.txt --gtf $PWDgtf/Homo_sapiens.Ensembl.GRCh37.72.gtf --bi /fdb/STAR_indices/2.7.8a/GENCODE/Gencode_human/release_27/genes-100 --od out_test1 -t paired --nthread $SLURM_CPUS_PER_TASK --readLength 50  --tophatAnchor 8 --cstat 0.0001 --tstat 6 --tmp /lscratch/${SLURM_JOB_ID}
rmats.py  --s1 $PWD s3.txt --s2 $PWD s4.txt --gtf $PWDgtf/Homo_sapiens.Ensembl.GRCh37.72.gtf --bi /fdb/STAR_indices/2.7.8a/GENCODE/Gencode_human/release_27/genes-100 --od out_test2 -t paired --nthread $SLURM_CPUS_PER_TASK --readLength 50  --tophatAnchor 8 --cstat 0.0001 --tstat 6 --tmp /lscratch/${SLURM_JOB_ID}
rmats.py  --s1 $PWD s5.txt --s2 $PWD s6.txt --gtf $PWDgtf/Homo_sapiens.Ensembl.GRCh37.72.gtf --bi /fdb/STAR_indices/2.7.8a/GENCODE/Gencode_human/release_27/genes-100 --od out_test3 -t paired --nthread $SLURM_CPUS_PER_TASK --readLength 50  --tophatAnchor 8 --cstat 0.0001 --tstat 6 --tmp /lscratch/${SLURM_JOB_ID}

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f rmats.swarm [-g 30] [-t 16] --module rmats
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module TEMPLATE Loads the TEMPLATE module for each subjob in the swarm 
 | |
 | |
 | |








