

document.querySelector('title').textContent = 'Plink on Biowulf';
Plink on Biowulf


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


Plink is a whole genome association analysis toolset, designed to perform a
range of basic, large-scale analyses in a computationally efficient manner.


The focus of PLINK is purely on analysis of genotype/phenotype data, so
there is no support for steps prior to this (e.g. study design and planning,
generating genotype calls from raw data). Through integration with [gPLINK](http://zzz.bwh.harvard.edu/plink/gplink.shtml) and
[Haploview](http://www.broad.mit.edu/mpg/haploview/), there is some
support for the subsequent visualization, annotation and storage of
results.


PLINK (one syllable) is being developed by Shaun Purcell at the Center for
Human Genetic Research ([CHGR](http://massgeneral.org/chgr/)),
Massachusetts General Hospital ([MGH](http://www.mgh.harvard.edu/)),
and the [Broad Institute](http://www.broad.mit.edu/) of Harvard
& MIT, with the [support of
others](http://zzz.bwh.harvard.edu/plink/contact.shtml#cite).



### References:


* Purcell et al., PLINK: A Tool Set for Whole-Genome Association and Population-Based Linkage Analyses. 2007. [Link](https://www.sciencedirect.com/science/article/pii/S0002929707613524)


Documentation
* [Plink Main Site](http://zzz.bwh.harvard.edu/plink/)


Important Notes
* Module Name: plink (see [the modules page](/apps/modules.html) for more information)
* Singlethreaded app
* Example files in /usr/local/apps/plink/TEST\_DATA



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

[user@cn3144 ~]$ **module load plink**
[+] Loading plink  3.6-alpha

[user@cn3144 ~]$ **cp /usr/local/apps/plink/TEST\_DATA/\* .**

[user@cn3144 ~]$ **plink2 --dummy 2 2 --freq --make-bed --out toy\_data**
PLINK v2.00a3.6LM 64-bit Intel (14 Aug 2022)   www.cog-genomics.org/plink/2.0/
(C) 2005-2022 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to toy_data.log.
Options in effect:
  --dummy 2 2
  --freq
  --make-bed
  --out toy_data

Start time: Wed Sep 21 09:35:39 2022
233597 MiB RAM detected; reserving 116798 MiB for main workspace.
Allocated 65698 MiB successfully, after larger attempt(s) failed.
Using up to 14 threads (change this with --threads).
Dummy data (2 samples, 2 SNPs) written to toy_data-temporary.pgen +
toy_data-temporary.pvar + toy_data-temporary.psam .
2 samples (2 females, 0 males; 2 founders) loaded from toy_data-temporary.psam.
2 variants loaded from toy_data-temporary.pvar.
1 binary phenotype loaded (2 cases, 0 controls).
Calculating allele frequencies... done.
--freq: Allele frequencies (founders only) written to toy_data.afreq .
Writing toy_data.fam ... done.
Writing toy_data.bim ... done.
Writing toy_data.bed ... done.
End time: Wed Sep 21 09:35:39 2022


[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. plink.sh). For example:



```

#!/bin/bash
cd /data/$USER/plink/t1
plink2  --dummy 2 2 --freq --make-bed --out toy_data

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --mem=6g plink.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. plink.swarm). For example:



```

cd /data/$USER/myseqs; plink2 --noweb --ped file1.ped --map file2.map --assoc
cd /data/$USER/myseqs; plink2 --noweb --ped file2.ped --map file2.map --assoc
cd /data/$USER/myseqs; plink2 --noweb --ped file3.ped --map file3.map --assoc
[...etc...]

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f plink.swarm [-g #] --module plink
```

where


|  |  |  |  |
| --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | --module plink Loads the plink module for each subjob in the swarm 
 | |
 | |








