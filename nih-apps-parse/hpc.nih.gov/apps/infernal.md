

document.querySelector('title').textContent = ' Infernal on Biowulf';
 Infernal on Biowulf


|  |
| --- |
| 
Quick Links
[Interactive Jobs](#int)
[Batch job on Biowulf](#batch)
[Using Swarm](#swarm)
[Documentation](#doc)
 |


Infernal is used to search sequence databases for homologs of structural RNA sequences, and to make
sequence- and structure-based RNA sequence alignments. Infernal builds a profile from a structurally
annotated multiple sequence alignment of an RNA family with a position-specific scoring system for substitutions, insertions, and deletions. Positions in the profile that are basepaired in the h consensus secondary
structure of the alignment are modeled as dependent on one another, allowing Infernal’s scoring system to
consider the secondary structure, in addition to the primary sequence, of the family being modeled. Infernal
profiles are probabilistic models called “covariance models”, a specialized type of stochastic context-free
grammar (SCFG) (Lari and Young, 1990).
Compared to other alignment and database search tools based only on sequence comparison, Infernal
aims to be significantly more accurate and more able to detect remote homologs because it models sequence and structure. But modeling structure comes at a high computational cost, and the slow speed of
CM homology searches has been a serious limitation of previous versions. With version 1.1, typical homology searches are now about 100x faster, thanks to the incorporation of accelerated HMM methods from the
HMMER3 software package (http://hmmer.org), making Infernal a much more practical tool for RNA
sequence analysis.


Documentation
<http://eddylab.org/infernal/>


Important Notes
* Module Name: epic2 (see [the modules page](/apps/modules.html) for more information)
* Example files are under /usr/local/apps/infernal/1.1.3/tutorial
* Multithreaded


Submitting an interactive job

Allocate an interactive session and run the interactive job there.



salloc.exe: Granted job allocation 789523
```
[biowulf]$ **sinteractive --mem=10g -cpus-per-task=4**
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn0135 are ready for job

[cn0135]$ **cd /data/$USER/**

[cn0135]$ **module load infernal**

[cn0135]$ **cat /usr/local/apps/infernal/1.1.3/tutorial/tRNA5.c.cm \
/usr/local/apps/infernal/1.1.3/tutorial/5S\_rRNA.c.cm \
/usr/local/apps/infernal/1.1.3/tutorial/Cobalamin.c.cm > minifam.cm**
[cn0135]$ **cmpress minifam.cm**
[cn0135]$ **cmscan --cpu 4 minifam.cm /usr/local/apps/infernal/1.1.3/tutorial/metag-example.fa**

[cn0135]$ **exit**
salloc.exe: Job allocation 789523 has been revoked.

[biowulf]$
```

The programs in infernal do not recognize $SLURM\_CPUS\_PER\_TASK variable.

Submitting a single batch job
1. Create a script file (myscript) similar to the one below

```
#! /bin/bash
# myscript
set -e

module load infernal || exit 1
cd /data/$USER/
cmscan --cpu 4 minifam.cm /usr/local/apps/infernal/1.1.3/tutorial/metag-example.fa

```

2. Submit the script on biowulf: 
 
 
```
[biowulf]$ sbatch --mem=10g --cpus-per-task=4 myscript
```


Using Swarm

Using the 'swarm' utility, one can submit many jobs to the cluster to run concurrently. 


Set up a swarm command file (eg /data/$USER/cmdfile). 



```

cd /data/$USER/dir1; cmscan ...
cd /data/$USER/dir2; cmscan ...
cd /data/$USER/dir3; cmscan ...
...
cd /data/$USER/dir20; cmscan ...

```


 submit the swarm job:
 
```
$ swarm -f cmdfile --module cmscan -g 10 -t 4
```

For more information regarding running swarm, see <swarm.html>





















