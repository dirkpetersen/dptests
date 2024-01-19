

document.querySelector('title').textContent = 'CRISPResso on Biowulf';
CRISPResso on Biowulf


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



CRISPResso2 is a software pipeline for the analysis of genome editing experiments. It is designed to enable rapid and intuitive interpretation of results produced by amplicon sequencing.
CRISPResso automatizes and performs the following steps summarized in the figure below:


1. filters low quality reads,
2. trims adapters,
3. aligns the reads to a reference amplicon,
4. quantifies the proportion of HDR and NHEJ outcomes,
5. quantifies frameshift/inframe mutations (if applicable) and identifies affected splice sites,
6. produces a graphical report to visualize and quantify the indels distribution and position.





### Reference:


* [Pinello, Luca, et al. "Analyzing CRISPR genome-editing experiments with CRISPResso." *Nature biotechnology* 34.7 (2016): 695-697.](https://www.nature.com/articles/nbt.3583)


Documentation
* [CRISPResso Main Site](https://crispresso.pinellolab.partners.org/submission)
* [Manual](https://crispresso.pinellolab.partners.org/help)
* [GitHub](https://github.com/pinellolab/CRISPResso2)


Important Notes
* CRISPResso and CRISPResso2 are both installed. Use module avail crispresso to see which versions are installed.
* Module Name: crispresso (see [the modules page](/apps/modules.html) for more information)



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. This example runs through the test suite supplied by the developer. Example commands can be found within the [testRelease.sh script](https://github.com/pinellolab/CRISPResso2/blob/master/tests/testRelease.sh).   
Sample session (user input in **bold**):



```

[user@biowulf crispresso]$ **sinteractive -c2 --mem=4g --gres=lscratch:10**
salloc.exe: Pending job allocation 11290667
salloc.exe: job 11290667 queued and waiting for resources
salloc.exe: job 11290667 has been allocated resources
salloc.exe: Granted job allocation 11290667
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn0863 are ready for job
srun: error: x11: no local DISPLAY defined, skipping
error: unable to open file /tmp/slurm-spank-x11.11290667.0
slurmstepd: error: x11: unable to read DISPLAY value

[user@cn0863 crispresso]$ **cd /lscratch/$SLURM\_JOB\_ID**

[user@cn0863 11290667]$ **git clone https://github.com/pinellolab/CRISPResso2.git**
Cloning into 'CRISPResso2'...
remote: Enumerating objects: 57, done.
remote: Counting objects: 100% (57/57), done.
remote: Compressing objects: 100% (44/44), done.
remote: Total 1118 (delta 29), reused 26 (delta 12), pack-reused 1061
Receiving objects: 100% (1118/1118), 1.53 MiB | 0 bytes/s, done.
Resolving deltas: 100% (783/783), done.

[user@cn0863 11290667]$ **cd CRISPResso2/tests/**

[user@cn0863 tests]$ **git checkout v2.2.14**
Note: checking out 'v2.2.14'.

You are in 'detached HEAD' state. You can look around, make experimental
changes and commit them, and you can discard any commits you make in this
state without impacting any branches by performing another checkout.

If you want to create a new branch to retain commits you create, you may
do so (now or later) by using -b with the checkout command again. Example:

  git checkout -b new_branch_name

HEAD is now at 4598226... HDR Updates - yw #82

[user@cn0863 tests]$ **module load crispresso**
[+] Loading crispresso  2.2.14  on cn0863
[+] Loading singularity  3.10.5  on cn0863

[user@cn0863 tests]$ **./testRelease.sh**
Running CRISPResso
Running CRISPResso with parameters
Running CRISPRessoBatch

[user@cn0863 tests]$ **exit**
exit
srun: error: cn0863: task 0: Exited with exit code 255
salloc.exe: Relinquishing job allocation 11290667
salloc.exe: Job allocation 11290667 has been revoked.

[user@biowulf crispresso]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. crispresso.sh). For example:



```

#!/bin/bash
set -e
module load crispresso
SEQ=CGGATGTTCCAATCAGTACGCAGAGAGTCGCCGTCTCCAAGGTGAAAGCGGAAGTAGGGCCTTCGCGCACCTCATGGAATCCCTTCTGCAGCACCTGGATCGCTTTTCCGAGCTTCTGGCGGTCTCAAGCACTACCTACGTCAGCACCTGGGACCCCGCCACCGTGCGCCGGGCCTTGCAGTGGGCGCGCTACCTGCGCCACATCCATCGGCGCTTTGGTCGG
CRISPResso -r1 FANC.Cas9.fastq -a $SEQ -g GGAATCCCTTCTGCAGCACC

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] crispresso.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. crispresso.swarm). For example:



```

CRISPResso -r1 A.fastq -a $SEQ -g GGAATCCCTTCTGCAGCACC
CRISPResso -r1 B.fastq -a $SEQ -g GGAATCCCTTCTGCAGCACC
CRISPResso -r1 C.fastq -a $SEQ -g GGAATCCCTTCTGCAGCACC
CRISPResso -r1 D.fastq -a $SEQ -g GGAATCCCTTCTGCAGCACC

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f crispresso.swarm [-g #] [-t #] --module crispresso
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module crispresso Loads the crispresso module for each subjob in the swarm 
 | |
 | |
 | |








