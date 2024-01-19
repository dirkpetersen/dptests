

document.querySelector('title').textContent = 'CAVIAR on Biowulf';
CAVIAR on Biowulf


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



CAVIAR (CAusal Variants Identication in Associated Regions), is a statistical framework that quantifies the probability of each variant to be causal while allowing with arbitrary number of causal variants.



### References:


* Farhad Hormozdiari, Ayellet V. Segre, Martijn van de Bunt, Xiao Li, Jong Wha J Joo, Michael Bilow, Jae Hoon Sul, Bogdan Pasaniuc and Eleazar Eski. Joint Fine Mapping of GWAS and eQTL Detects Target Gene and Relevant Tissue. The American Journal of Human Genetics (AJHG). [doi:10.1016/j.ajhg.2016.10.003](https://dx.doi.org/10.1016/j.ajhg.2016.10.003)
* Farhad Hormozdiari, Emrah Kostem, Eun Yong Kang, Bogdan Pasaniuc and Eleazar Eskin. Identifying Causal Variants at Loci with Multiple Signals of Association. Genetics, 44, 725–731 (2014). [doi:10.1534/genetics.114.167908](https://doi.org/10.1534/genetics.114.167908)


Documentation
* [CAVIAR Main Site](http://genetics.cs.ucla.edu/caviar/index.html)


Important Notes
* Module Name: caviar (see [the modules page](/apps/modules.html) for more information)
* Multithreaded/singlethreaded/MPI...
* environment variables set
	+ CAVIAR\_HOME* Example files in $CAVIAR\_HOME/sample\_data



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --gres lscratch:1**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job
[user@cn3144 ~]$ **module load caviar**
[+] Loading caviar, version a97e614...
[user@cn3144 ~]$ **cd /lscratch/$SLURM\_JOBID**
[user@cn3144 46116226]$ **cp -r $CAVIAR\_HOME/sample\_data/\* .**
[user@cn3144 46116226]$ **CAVIAR -l DDB1.top100.sig.SNPs.ld -z DDB1.top100.sig.SNPs.ZScores -o outfile**
@-------------------------------------------------------------@
| CAVIAR!|    v2.0         |  01/Aug/2017 | 
|-------------------------------------------------------------|
| (C) 2017 Farhad Hormozdiari, GNU General Public License, v2 |
|-------------------------------------------------------------|
| For documentation, citation & bug-report instructions:      |
| http://genetics.cs.ucla.edu/caviar/            |
@-------------------------------------------------------------@
00
add
1.73212e-750.1
FINISH
reach=100
100
Max Causal=2
Total=5051
0 1.36408e-06
Total Likelihood= 4.552733e+48 SNP=100 

56 5.000000e-01
53 9.999984e-01

[user@cn3144 46116226]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. caviar.sh). For example:



```

#!/bin/sh

set -e

module load caviar
cp -r $CAVIAR_HOME/sample_data/* .
CAVIAR -l DDB1.top100.sig.SNPs.ld -z DDB1.top100.sig.SNPs.ZScores -o outfile

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] caviar.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. caviar.swarm). For example:



```

CAVIAR -l ldfile1 -z zfile1 -o outfile1
CAVIAR -l ldfile2 -z zfile2 -o outfile2
CAVIAR -l ldfile3 -z zfile3 -o outfile3
CAVIAR -l ldfile4 -z zfile4 -o outfile4

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f caviar.swarm [-g #] [-t #] --module caviar
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g #  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t #  Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module caviar  Loads the TEMPLATE module for each subjob in the swarm 
 | |
 | |
 | |








