

document.querySelector('title').textContent = 'Flashpca on Biowulf';
Flashpca on Biowulf


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



FlashPCA performs fast principal component analysis (PCA) of single nucleotide polymorphism (SNP) data, similar to smartpca from EIGENSOFT (http://www.hsph.harvard.edu/alkes-price/software/) and shellfish (https://github.com/dandavison/shellfish). FlashPCA is based on the https://github.com/yixuan/spectra/ library.



Documentation
* [Flashpca Main Site](https://github.com/gabraham/flashpca)


Important Notes
* Module Name: flashpca (see [the modules page](/apps/modules.html) for more information)
* singlethreaded app
* Example files in /usr/local/apps/flashpca/TESTDATA/



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load flashpca**

[user@cn3144 ~]$ **cp /usr/local/apps/flashpca/TESTDATA/merged\* .**

[user@cn3144 ~]$ **flashpca --bed merged.bed --fam merged.fam --bim merged.bim**
[Mon Sep 24 11:57:53 2018] Start flashpca (version 2.0)
[Mon Sep 24 11:57:53 2018] seed: 1
[Mon Sep 24 11:57:53 2018] blocksize: 24766 (264699008 bytes per block)
[Mon Sep 24 11:57:53 2018] PCA begin
[Mon Sep 24 11:58:00 2018] PCA done
[Mon Sep 24 11:58:00 2018] Writing 10 eigenvalues to file eigenvalues.txt
[Mon Sep 24 11:58:00 2018] Writing 10 eigenvectors to file eigenvectors.txt
[Mon Sep 24 11:58:00 2018] Writing 10 PCs to file pcs.txt
[Mon Sep 24 11:58:00 2018] Writing 10 proportion variance explained to file pve.txt
[Mon Sep 24 11:58:00 2018] Goodbye!

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. flashpca.sh). For example:



```

#!/bin/bash
set -e
module load flashpca
flashpca --bed merged.bed --fam merged.fam --bim merged.bim

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch flashpca.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. flashpca.swarm). For example:



```

cd folder1; flashpca --bed merged.bed --fam merged.fam --bim merged.bim
cd folder2; flashpca --bed merged.bed --fam merged.fam --bim merged.bim
cd folder3; flashpca --bed merged.bed --fam merged.fam --bim merged.bim

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f flashpca.swarm --module flashpca
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module flashpca Loads the flashpca module for each subjob in the swarm 
 | |
 | |
 | |








