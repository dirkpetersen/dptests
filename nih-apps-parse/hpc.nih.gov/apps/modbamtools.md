

document.querySelector('title').textContent = 'modbamtools on Biowulf';
modbamtools on Biowulf


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



Modbamtools is a set of tools to manipulate and visualize DNA/RNA base modification data that are stored in bam format. htslib has included a support for parsing modified base tags from alignment files (MM and ML). These tags have provided a better/efficient way for storing modification data inside alignment files. 








### References:


* ProfileRoham Razaghi, ProfilePaul W. Hook, ProfileShujun Ou, ProfileMichael C. Schatz, ProfileKasper D. Hansen, ProfileMiten Jain, ProfileWinston Timp
*Modbamtools: Analysis of single-molecule epigenetic data for long-range profiling, heterogeneity, and clustering*
[Biorxiv](https://www.biorxiv.org/content/10.1101/2022.07.07.499188v1.article-info)


Documentation
* modbamtools Main Site:[Main Site](https://rrazaghi.github.io/modbamtools/)
* modbamtools Main Site:[Tutorial](https://rrazaghi.github.io/modbamtools/tutorial/)


Important Notes
* Module Name: modbamtools (see [the modules page](/apps/modules.html) for more information)
* Environment variables set 
	+ $MODBAMTOOLS\_TEST\_DATA #include bam, bigwig, gtf and bed files.



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive --cpus-per-task=2 --mem=2G**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ]$ **module load modbamtools**
[user@cn3144 ]$ **cp -r ${MODBAMTOOLS\_TEST\_DATA:-none}/\* .**
[user@cn3144 ]$ **cd modbamtools\_tutorial\_files**
[user@cn3144 ]$ **modbamtools plot -r chr20:58815000-58895000 \
 --gtf gencode.v38.annotation.sorted.gtf.gz \
 --out . \
 --prefix gm12878\_GNAS \
 --samples GM12878 \
 --track-titles Genes\
 gm12878\_ul\_sup\_megalodon\_HP\_chr20.bam**
[user@cn3144]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```

This will generate gm12878\_GNAS.html:


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. modbamtools.sh). For example:



```

#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G
#SBATCH --time=2:00:00
#SBATCH --partition=norm

set -e
module load modbamtools
cp -r ${MODBAMTOOLS_TEST_DATA:-none}/* .
cd modbamtools_tutorial_files
modbamtools plot -r chr20:58815000-58895000 \
    --gtf gencode.v38.annotation.sorted.gtf.gz \
    --out . \
    --hap \
    --prefix gm12878_GNAS \
    --samples GM12878 \
    --track-titles Genes\
    gm12878_ul_sup_megalodon_HP_chr20.bam 

```

 Submit the job:

```
sbatch modbamtools.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. job.swarm). For example:



```

       modbamtools plot xxx
       modbamtools plot xxx
    
```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f job.swarm [-g #] --module modbamtools
```

where
 

|  |  |
| --- | --- |
| -g *#*  | Number of Gigabytes of memory required for each process
 (1 line in the swarm command file)  |
| --module  | Loads the module for each subjob in the swarm  |










