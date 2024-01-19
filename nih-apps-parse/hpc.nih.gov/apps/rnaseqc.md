

document.querySelector('title').textContent = ' RNA-SeQC on Biowulf';
 RNA-SeQC on Biowulf


|  |
| --- |
| 
Quick Links
[Interactive Jobs](#int)
[Batch job on Biowulf](#batch)
[Using Swarm](#swarm)
[Documentation](#doc)
 |


Fast, efficient RNA-Seq metrics for quality control and process optimization.

Documentation
<https://github.com/getzlab/rnaseqc>


Important Notes
* Module Name: rnaseqc (see [the modules page](/apps/modules.html) for more information)
* Example files are under /usr/local/apps/rnaseqc/test\_data


Submitting an interactive job

Allocate an interactive session and run the interactive job there.



```
[biowulf]$ **sinteractive --mem=5g**
salloc.exe: Granted job allocation 789523
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn0135 are ready for job

[cn0135]$ **module load rnaseqc**

[cn0135]$ **cp -r /usr/local/apps/rnaseqc/test\_data /data/$USER/**
[cn0135]$ **cd /data/$USER/test\_data**
[cn0135]$ **rnaseqc chr1.gtf chr1.bam user.out**
[cn0135]$ **exit**
salloc.exe: Job allocation 789523 has been revoked.

[biowulf]$
```


Submitting a single batch job
1. Create a script file (myscript) similar to the one below

```
#! /bin/bash
# myscript
set -e

module load rnaseqc || exit 1
cd /data/$USER/test_data
rnaseqc chr1.gtf chr1.bam user.out

```

2. Submit the script on biowulf: 
 
 
```
[biowulf]$ sbatch --mem=5g myscript
```


Using Swarm

Using the 'swarm' utility, one can submit many jobs to the cluster to run concurrently. 


Set up a swarm command file (eg /data/$USER/cmdfile). 



```

cd /data/$USER/dir1; rnaseqc chr1.gtf chr1.bam user.out
cd /data/$USER/dir2; rnaseqc chr1.gtf chr1.bam user.out
cd /data/$USER/dir3; rnaseqc chr1.gtf chr1.bam user.out
...
cd /data/$USER/dir20; rnaseqc chr1.gtf chr1.bam user.out

```


 submit the swarm job:
 
```
$ swarm -f cmdfile --module rnaseqc -g 5
```

For more information regarding running swarm, see <swarm.html>



















