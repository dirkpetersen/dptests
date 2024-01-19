

document.querySelector('title').textContent = ' Epic2 on Biowulf';
 Epic2 on Biowulf


|  |
| --- |
| 
Quick Links
[Interactive Jobs](#int)
[Batch job on Biowulf](#batch)
[Using Swarm](#swarm)
[Documentation](#doc)
 |


epic2 is an ultraperformant reimplementation of SICER. It focuses on speed, low memory overhead and ease of use.

It also contains a reimplementation of the SICER-df scripts for differential enrichment and a script to create many kinds of bigwigs from your data.
**Features:**
* easy to install and use
* reads sam, single-end bam, bed and bedpe (.gz)
* extremely fast
* very low memory requirements
* works both with and without input
* metadata for ~80 UCSC genomes built in
* easily use custom genomes and assemblies with --chromsizes and --effective-genome-fraction args
* differential enrichment for WT vs. KO works with and without input (epic2-df)
* fixes two bugs in the original SICER and one bug in epic
* create many types of useful bigwigs for visualization in genome browsers


Documentation
<https://github.com/biocore-ntnu/epic2>


Important Notes
* Module Name: epic2 (see [the modules page](/apps/modules.html) for more information)
* Example files are under /usr/local/apps/epic2/examples


Submitting an interactive job

Allocate an interactive session and run the interactive job there.



```
[biowulf]$ **sinteractive --mem=5g**
salloc.exe: Granted job allocation 789523
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn0135 are ready for job

[cn0135]$ **cd /data/$USER/**

[cn0135]$ **module load epic2**

[cn0135]$ **epic2 -ex**
# An example of command can be copied from the output.
[cn0135]$ **epic2 -t /opt/conda/envs/app/lib/python3.7/site-packages/epic2/examples/test.bed.gz \
-c /opt/conda/envs/app/lib/python3.7/site-packages/epic2/examples/control.bed.gz > deleteme.txt**

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

module load epic2 || exit 1
cd /data/$USER/
epic2 -t /opt/conda/envs/app/lib/python3.7/site-packages/epic2/examples/test.bed.gz \
-c /opt/conda/envs/app/lib/python3.7/site-packages/epic2/examples/control.bed.gz > deleteme.txt
```

2. Submit the script on biowulf: 
 
 
```
[biowulf]$ sbatch --mem=5g myscript
```


Using Swarm

Using the 'swarm' utility, one can submit many jobs to the cluster to run concurrently. 


Set up a swarm command file (eg /data/$USER/cmdfile). 



```

cd /data/$USER/dir1; epic2 ...
cd /data/$USER/dir2; epic2 ...
cd /data/$USER/dir3; epic2 ...
...
cd /data/$USER/dir20; epic2 ...

```


 submit the swarm job:
 
```
$ swarm -f cmdfile --module epic2 -g 5
```

For more information regarding running swarm, see <swarm.html>





















