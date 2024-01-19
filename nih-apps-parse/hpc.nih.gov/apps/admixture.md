

document.querySelector('title').textContent = 'Admixture on Biowulf';
Admixture on Biowulf


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



ADMIXTURE is a software tool for maximum likelihood estimation of individual ancestries from multilocus SNP genotype datasets. It uses the same statistical model as STRUCTURE but calculates estimates much more rapidly using a fast numerical optimization algorithm. 


### References:


* D.H. Alexander, J. Novembre, and K. Lange. Fast model-based estimation of ancestry in unrelated individuals. Genome Research, 19:1655–1664, 2009.


Documentation
* [Admixture Main Site](https://dalexander.github.io/admixture/index.html)


Important Notes
* Module Name: admixture (see [the modules page](/apps/modules.html) for more information)
* Multithreaded app (use -jX option)
* Example files in /usr/local/apps/admixture/example



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session (user input in **bold**):



```

[user@biowulf ~]$ **sinteractive -c4 --mem=8g --gres=lscratch:10**
salloc.exe: Pending job allocation 11085118
salloc.exe: job 11085118 queued and waiting for resources
salloc.exe: job 11085118 has been allocated resources
salloc.exe: Granted job allocation 11085118
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn0848 are ready for job
srun: error: x11: no local DISPLAY defined, skipping
error: unable to open file /tmp/slurm-spank-x11.11085118.0
slurmstepd: error: x11: unable to read DISPLAY value

[user@cn0848 ~]$ **cd /lscratch/$SLURM\_JOB\_ID**

[user@cn0848 11085118]$ **cp /usr/local/apps/admixture/example/\* .**

[user@cn0848 11085118]$ **module load admixture**
[+] Loading admixture  1.3.0  on cn0848

[user@cn0848 11085118]$ **admixture -help**
****                   ADMIXTURE Version 1.3.0                  ****
****                    Copyright 2008-2015                     ****
****           David Alexander, Suyash Shringarpure,            ****
****                John  Novembre, Ken Lange                   ****
****                                                            ****
****                 Please cite our paper!                     ****
****   Information at www.genetics.ucla.edu/software/admixture  ****

Block size set manually to 0 SNPs
Usage: admixture <input file> <K>
See --help or manual for more advanced usage.

[user@cn0848 11085118]$ **admixture hapmap3.bed 3 -j${SLURM\_CPUS\_PER\_TASK}**
****                   ADMIXTURE Version 1.3.0                  ****
****                    Copyright 2008-2015                     ****
****           David Alexander, Suyash Shringarpure,            ****
****                John  Novembre, Ken Lange                   ****
****                                                            ****
****                 Please cite our paper!                     ****
****   Information at www.genetics.ucla.edu/software/admixture  ****

Random seed: 43
Point estimation method: Block relaxation algorithm
Convergence acceleration algorithm: QuasiNewton, 3 secant conditions
Point estimation will terminate when objective function delta < 0.0001
Estimation of standard errors disabled; will compute point estimates only.
Size of G: 324x13928
Performing five EM steps to prime main algorithm
1 (EM)  Elapsed: 0.3    Loglikelihood: -4.38757e+06     (delta): 2.87325e+06
2 (EM)  Elapsed: 0.3    Loglikelihood: -4.25681e+06     (delta): 130762
3 (EM)  Elapsed: 0.299  Loglikelihood: -4.21622e+06     (delta): 40582.9
4 (EM)  Elapsed: 0.299  Loglikelihood: -4.19347e+06     (delta): 22748.2
5 (EM)  Elapsed: 0.299  Loglikelihood: -4.17881e+06     (delta): 14663.1
Initial loglikelihood: -4.17881e+06
Starting main algorithm
1 (QN/Block)    Elapsed: 0.723  Loglikelihood: -3.94775e+06     (delta): 231058
[...snip]
Summary:
Converged in 21 iterations (21.566 sec)
Loglikelihood: -3799887.171935
Fst divergences between estimated populations:
        Pop0    Pop1
Pop0
Pop1    0.163
Pop2    0.073   0.156
Writing output files.

[user@cn0848 11085118]$ **exit**
exit
salloc.exe: Relinquishing job allocation 11085118

[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. admixture.sh). For example:



```

#!/bin/bash
set -e
module load admixture
admixture hapmap3.bed 3 -j${SLURM_CPUS_PER_TASK} > admixture.out

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=16 --mem=20g admixture.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. admixture.swarm). For example:



```

admixture hapmap3.bed 3 -j${SLURM_CPUS_PER_TASK} > admixture3.out
admixture hapmap20.bed 20 -j${SLURM_CPUS_PER_TASK} > admixture20.out
admixture hapmap1000.bed 1000 -j${SLURM_CPUS_PER_TASK} > admixture1000.out


```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f admixture.swarm -g 20 -t 16 --module admixture
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module admixture Loads the admixture module for each subjob in the swarm 
 | |
 | |
 | |








