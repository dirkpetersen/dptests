

document.querySelector('title').textContent = 'PBSuite on Biowulf';
PBSuite on Biowulf


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



The PBSuite contains two projects created for analysis of Pacific Biosciences long-read sequencing data: PBHoney and PBJelly.




PBHoney is an implementation of two variant-identification approaches designed to exploit the high mappability of long reads (i.e., greater than 10,000 bp). PBHoney considers both intra-read discordance and soft-clipped tails of long reads to identify structural variants.




PBJelly is a highly automated pipeline that aligns long sequencing reads (such as PacBio RS reads or long 454 reads in fasta format) to high-confidence draft assembles. PBJelly fills or reduces as many captured gaps as possible to produce upgraded draft genomes.



### References:


* English AC, Richards S, Han Y, Wang M, Vee V, et al. (2012) Mind the Gap: Upgrading Genomes with Pacific Biosciences RS Long-Read Sequencing Technology. PLOS ONE 7(11): e47768. [doi:10.1371/journal.pone.0047768](https://doi.org/10.1371/journal.pone.0047768)
* English A.C., Salerno W.J. & Reid J.G. (2014) PBHoney: identifying genomic variants via long-read discordance and interrupted mapping. BMC Bioinformatics 15, 180. [doi:10.1186/1471-2105-15-180](https://doi.org/10.1186/1471-2105-15-180)


Documentation
* [PBSuite Homepage](https://sourceforge.net/projects/pb-jelly/)
* $PBSUITE\_HOME/docs


Important Notes
* Module Name: pbsuite or pbjelly (see [the modules page](/apps/modules.html) for more information)
* Environment variables set 
	+ PBSUITE\_HOME* Example files in $PBSUITE\_HOME/docs/honeyExample/ and $PBSUITE\_HOME/docs/jellyExample/



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample sessions (user input in **bold**):


### PBHoney



```

[user@biowulf]$ **sinteractive**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load pbsuite**
[user@cn3144 ~]$ **cp $PBSUITE\_HOME/docs/honeyExample/\* .**
[user@cn3144 ~]$ **sh workflow.sh**
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```

### PBJelly




```

[user@biowulf]$ **sinteractive**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load pbsuite**
[user@cn3144 ~]$ **cp -r $PBSUITE\_HOME/docs/jellyExample/\* .**
[user@cn3144 ~]$ **sed -i "s|/\_\_PATH\_\_/\_TO\_/jellyExample|$PWD|g" Protocol.xml**
[user@cn3144 ~]$ **for stage in setup mapping support extraction assembly output
do
 Jelly.py $stage Protocol.xml
done**
[user@cn3144 ~]$ **summarizeAssembly.py jelly.out.fasta** *# check your results*
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. pbjelly.sh). For example:



```

#!/bin/sh
set -e
module load pbsuite

for stage in setup mapping support extraction assembly output
do
    Jelly.py $stage Protocol.xml
done

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] pbjelly.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. pbjelly.swarm). For example:



```

for stage in setup mapping support extraction assembly output; do Jelly.py $stage protocol1.xml || exit 1; done
for stage in setup mapping support extraction assembly output; do Jelly.py $stage protocol2.xml || exit 1; done
for stage in setup mapping support extraction assembly output; do Jelly.py $stage protocol3.xml || exit 1; done
for stage in setup mapping support extraction assembly output; do Jelly.py $stage protocol4.xml || exit 1; done

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f pbjelly.swarm [-g #] [-t #] --module pbsuite
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module pbsuite Loads the pbsuite module for each subjob in the swarm 
 | |
 | |
 | |










