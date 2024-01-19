

document.querySelector('title').textContent = ' Sickle on Biowulf';
 Sickle on Biowulf


|  |
| --- |
| 
Quick Links
[Interactive Jobs](#int)
[Batch job on Biowulf](#batch)
[Using Swarm](#swarm)
 |


Sickle is a windowed adaptive trimming tool for FASTQ files using quality. Most modern sequencing technologies produce reads that have deteriorating quality towards the 3'-end and some towards the 5'-end as well. Incorrectly called bases in both regions negatively impact assembles, mapping, and downstream bioinformatics analyses. Sickle is a tool that uses sliding windows along with quality and length thresholds to determine when quality is sufficiently low to trim the 3'-end of reads and also determines when the quality is sufficiently high enough to trim the 5'-end of reads. It will also discard reads based upon the length threshold. It takes the quality values and slides a window across them whose length is 0.1 times the length of the read. If this length is less than 1, then the window is set to be equal to the length of the read. Otherwise, the window slides along the quality values until the average quality in the window rises above the threshold, at which point the algorithm determines where within the window the rise occurs and cuts the read and quality there for the 5'-end cut. Then when the average quality in the window drops below the threshold, the algorithm determines where in the window the drop occurs and cuts both the read and quality strings there for the 3'-end cut. However, if the length of the remaining sequence is less than the minimum length threshold, then the read is discarded entirely. 5'-end trimming can be disabled.



Documentation
<https://github.com/najoshi/sickle/>


Important Notes
* Module Name: sickle (see [the modules page](/apps/modules.html) for more information)
* Example files are under /usr/local/apps/sickle/1.33/TEST\_DATA or $SICKLE\_TEST\_DATA


Submitting an interactive job

Allocate an interactive session and run the interactive job there.



```

[biowulf]$ **sinteractive --mem=5g**
salloc.exe: Granted job allocation 789523
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn0135 are ready for job

[cn0135]$ **cd /data/$USER/**

[cn0135]$ **module load sickle**

[cn0135]$ **cp -r $SICKLE\_TEST\_DATA .**

[cn0135]$ **cd TEST\_DATA**

[cn0135]$ **sickle pe -f test.f.fastq -r test.r.fastq -t sanger -o output1.fastq -p output2.fastq -s trimmed\_singles\_file.fastq**

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

module load sickle || exit 1
cd /data/$USER/test/
sickle pe -f test.f.fastq \
          -r test.r.fastq \
          -t sanger \
          -o output1.fastq \
          -p output2.fastq \
          -s trimmed_singles_file.fastq

```

2. Submit the script on biowulf: 
 
 
```
[biowulf]$ sbatch --mem=5g myscript
```


Using Swarm

Using the 'swarm' utility, one can submit many jobs to the cluster to run concurrently. 


Set up a swarm command file (eg /data/$USER/cmdfile). 



```

cd /data/$USER/dir1; sickle ...
cd /data/$USER/dir2; sickle ...
cd /data/$USER/dir3; sickle ...
...
cd /data/$USER/dir20; sickle ...

```


 submit the swarm job:
 
```
$ swarm -f cmdfile --module sickle -g 5
```

For more information regarding running swarm, see <swarm.html>



















