

document.querySelector('title').textContent = ' Augustus on Biowulf';
 Augustus on Biowulf


|  |
| --- |
| 
Quick Links
[Interactive Jobs](#int)
[Batch job on Biowulf](#batch)
[Using Swarm](#swarm)
[Documentation](#doc)
 |



AUGUSTUS is a gene prediction program written by Mario Stanke, Oliver Keller, Stefanie König, Lizzy Gerischer and Katharina Hoff. It can be used as an ab initio program, which means it bases its prediction purely on the sequence. AUGUSTUS may also incorporate hints on the gene structure coming from extrinsic sources such as EST, MS/MS, protein alignments and syntenic genomic alignments. Since version 3.0 AUGUSTUS can also predict the genes simultaneously in several aligned genomes (see README-cgp.txt). 

Documentation
<https://github.com/Gaius-Augustus/Augustus>


Important Notes
* Module Name: augustus (see [the modules page](/apps/modules.html) for more information)
* Example files are under /usr/local/apps/augustus/version/example\*


Submitting an interactive job

Allocate an interactive session and run the interactive job there.



```
[biowulf]$ **sinteractive --mem=10g**
salloc.exe: Granted job allocation 789523
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn0135 are ready for job

[cn0135]$ **cd /data/$USER/**

[cn0135]$ **module load augustus**

[cn0135]$ **cp -rp /usr/local/apps/augustus/version/examples/example.fa /data/$USER**
[cn0135]$ **augustus --species=human example.fa**

[cn0135]$ **exit**
salloc.exe: Job allocation 789523 has been revoked.

```


Submitting a single batch job
1. Create a script file (myscript) similar to the one below

```
#! /bin/bash
# myscript
set -e

module load augustus || exit 1
cd /data/$USER/
augustus --species=human example.fa
```

2. Submit the script on biowulf: 
 
 
```
[biowulf]$ sbatch --mem=5g myscript
```


Using Swarm

Using the 'swarm' utility, one can submit many jobs to the cluster to run concurrently. 


Set up a swarm command file (eg /data/$USER/cmdfile). 



```

cd /data/$USER/dir1; augustus species=human example.fa
cd /data/$USER/dir2; augustus species=human example.fa
cd /data/$USER/dir3; augustus species=human example.fa
...
cd /data/$USER/dir20; augustus species=human example.fa

```


 submit the swarm job:
 
```
$ swarm -f cmdfile --module augustus -g 5
```

For more information regarding running swarm, see <swarm.html>

















