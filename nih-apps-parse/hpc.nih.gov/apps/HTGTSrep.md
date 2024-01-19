

document.querySelector('title').textContent = 'HTGTSrep on Biowulf';
HTGTSrep on Biowulf


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



A pipeline for comprehensive analysis of HTGTS-Rep-seq.



Documentation
* [HTGTSrep on GitHub](https://github.com/Yyx2626/HTGTSrep)
* [HTGTSrep on BitBucket](https://bitbucket.org/adugduzhou/htgtsrep/src/master/)


Important Notes
* Module Name: HTGTSrep (see [the modules page](/apps/modules.html) for more information)
 * This application is run on Biowulf using a staff-written wrapper script. Simply run the command htgtsrep after loading the module. Do not try to invoke the HTGTSrep.py script using python.



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf ~]$ **sinteractive --cpus-per-task=4 --mem=8g --gres=lscratch:10**
salloc.exe: Pending job allocation 51996182
salloc.exe: job 51996182 queued and waiting for resources
salloc.exe: job 51996182 has been allocated resources
salloc.exe: Granted job allocation 51996182
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn0861 are ready for job
srun: error: x11: no local DISPLAY defined, skipping

[user@cn0861 ~]$ **git clone https://github.com/Yyx2626/HTGTSrep.git /lscratch/$SLURM\_JOB\_ID**
Cloning into '/lscratch/51996182'...
remote: Enumerating objects: 4, done.
remote: Counting objects: 100% (4/4), done.
remote: Compressing objects: 100% (4/4), done.
remote: Total 1203 (delta 0), reused 0 (delta 0), pack-reused 1199
Receiving objects: 100% (1203/1203), 485.90 MiB | 51.51 MiB/s, done.
Resolving deltas: 100% (324/324), done.
Checking out files: 100% (1226/1226), done.

[user@cn0861 ~]$ **cd /lscratch/$SLURM\_JOB\_ID/HTGTSrep/test**

[user@cn0861 test]$ **module load HTGTSrep**
[+] Loading HTGTSrep  9fe74ff  on cn0861
[+] Loading singularity  3.5.3  on cn0861

[user@cn0861 test]$ **htgtsrep run -m metadata.txt -r1 r1.fq.gz -r2 r2.fq.gz -o out**
Welcome to HTGTSrep pipeline!!!
INFO  @ Fri, 13 Mar 2020 16:40:31: Parameters: /opt/HTGTSrep/HTGTSrep/HTGTSrep.py run -m metadata.txt -r1 r1.fq.gz -r2 r2.fq.gz -o out
INFO  @ Fri, 13 Mar 2020 16:40:31: Parsing meta files.....
INFO  @ Fri, 13 Mar 2020 16:40:31: Preprocessing read files.....
INFO  @ Fri, 13 Mar 2020 16:40:32: fastq-multx -m 0 -x -b -d 0 -B out/barcodes.txt out/raw.r1.fq out/raw.r2.fq -o out/%_R1.fq out/%_R2.fq
Using Barcode File: out/barcodes.txt
Id      Count   File(s)
HC091_Alt231    16333   out/HC091_Alt231_R1.fq  out/HC091_Alt231_R2.fq
HC092_Alt231    13145   out/HC092_Alt231_R1.fq  out/HC092_Alt231_R2.fq
unmatched       0       out/unmatched_R1.fq     out/unmatched_R2.fq
total   29478
INFO  @ Fri, 13 Mar 2020 16:40:32: Joining reads......
[snip...]

[user@cn0861 test]$ **ls out/**
barcodes.txt  HC091_Alt231  HC092_Alt231  logs  metadata.txt  stat

[user@cn0861 test]$ **cp -r out/ /data/$USER/**

[user@cn0861 test]$ **exit**
exit
salloc.exe: Relinquishing job allocation 51996182
salloc.exe: Job allocation 51996182 has been revoked.

[user@biowulf ~]$ 

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. HTGTSrep.sh). For example:



```

#!/bin/bash
set -e
module load HTGTSrep
cd /data/$USER/HTGTSrep/test
htgtsrep run -m metadata.txt -r1 r1.fq.gz -r2 r2.fq.gz -o /data/$USER/out

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] HTGTSrep.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. HTGTSrep.swarm). For example:



```

htgtsrep run -m metadata.txt -r1 r1A.fq.gz -r2 r2A.fq.gz -o /data/$USER/outA
htgtsrep run -m metadata.txt -r1 r1B.fq.gz -r2 r2B.fq.gz -o /data/$USER/outB
htgtsrep run -m metadata.txt -r1 r1C.fq.gz -r2 r2C.fq.gz -o /data/$USER/outC
htgtsrep run -m metadata.txt -r1 r1D.fq.gz -r2 r2D.fq.gz -o /data/$USER/outD

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f HTGTSrep.swarm [-g #] [-t #] --module HTGTSrep
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module HTGTSrep Loads the HTGTSrep module for each subjob in the swarm 
 | |
 | |
 | |








