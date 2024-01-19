

document.querySelector('title').textContent = 'Trimgalore on Biowulf';
Trimgalore on Biowulf


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


Trim Galore! is a wrapper script to automate quality and adapter trimming as well as quality control, with some added functionality to remove biased methylation positions for RRBS sequence files (for directional, non-directional (or paired-end) sequencing). It's main features are:


* For adapter trimming, Trim Galore! uses the first 13 bp of Illumina standard adapters ('AGATCGGAAGAGC') by default (suitable for both ends of paired-end libraries), but accepts other adapter sequence, too
* For MspI-digested RRBS libraries, Trim Galore! performs quality and adapter trimming in two subsequent steps. This allows it to remove 2 additional bases that contain a cytosine which was artificially introduced in the end-repair step during the library preparation
* For any kind of FastQ file other than MspI-digested RRBS, Trim Galore! can perform single-pass adapter- and quality trimming
* The Phred quality of basecalls and the stringency for adapter removal can be specified individually
* Trim Galore! can remove sequences if they become too short during the trimming process. For paired-end files Trim Galore! removes entire sequence pairs if one (or both) of the two reads became shorter than the set length cutoff. Reads of a read-pair that are longer than a given threshold but for which the partner read has become too short can optionally be written out to single-end files. This ensures that the information of a read pair is not lost entirely if only one read is of good quality
* Trim Galore! can trim paired-end files by 1 additional bp from the 3' end of all reads to avoid problems with invalid alignments with Bowtie 1
* Trim Galore! accepts and produces standard or gzip compressed FastQ files
* FastQC can be run on the resulting output files once trimming has completed (optional)

 
Documentation * <http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/>



Important Notes * Module Name: trimgalore (see [the 
 modules page](/apps/modules.html) for more information)
* Example files in /usr/local/apps/trimgalore/version/test\_files





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

[user@cn3144 ~]$ **module load trimgalore**
[user@cn3144 ~]$ **trim\_galore inputfile**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```




Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. trimgalore.sh). For example:



```

#!/bin/bash
set -e
module load trimgalore
trim_galore inputfile [options...]
```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--mem=#] trimgalore.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. trimgalore.swarm). For example:



```

cd dir1; trim_galore infile
cd dir2; trim_galore infile
cd dir3; trim_galore infile
cd dir4; trim_galore infile

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f TEMPLATE.swarm [-g #] --module trimgalore
```

where
 

|  |  |
| --- | --- |
| -g *#*  | Number of Gigabytes of memory required for each process 
 (1 line in the swarm command file)  |
| --module  | Loads the module for each subjob in the swarm  |




