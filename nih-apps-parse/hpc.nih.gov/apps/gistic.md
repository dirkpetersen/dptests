

document.querySelector('title').textContent = 'GISTIC on Biowulf';
GISTIC on Biowulf


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



GISTIC, or Genomic Identification of Significant Targets in Cancer, identifies regions of the genome that are significantly amplified or deleted across a set of samples. Each aberration is assigned a G-score that considers the amplitude of the aberration as well as the frequency of its occurrence across samples. False Discovery Rate q-values are then calculated for the aberrant regions, and regions with q-values below a user-defined threshold are considered significant. For each significant region, a “peak region” is identified, which is the part of the aberrant region with greatest amplitude and frequency of alteration. In addition, a “wide peak” is determined using a leave-one-out algorithm to allow for errors in the boundaries in a single sample. The “wide peak” boundaries are more robust for identifying the most likely gene targets in the region. Each significantly aberrant region is also tested to determine whether it results primarily from broad events (longer than half a chromosome arm), focal events, or significant levels of both. The GISTIC module reports the genomic locations and calculated q-values for the aberrant regions. It identifies the samples that exhibit each significant amplification or deletion, and it lists genes found in each “wide peak” region.



### References:


* Mermel C, Schumacher S, et al. (2011). "GISTIC2.0 facilitates sensitive and confident localization of the targets of focal somatic copy-number alteration in human cancers." Genome Biology, 12:R41. [doi:10.1186/gb-2011-12-4-r41](https://doi.org/10.1186/gb-2011-12-4-r41)
* Beroukhim R, Mermel C, et al. (2010). "The landscape of somatic copy -number alteration across human cancers." Nature, 463:899-905. [doi:10.1038/nature08822](https://doi.org/10.1038/nature08822)
* Beroukhim R, Getz G, et al. (2007). “Assessing the significance of chromosomal abberations in cancer: Methodology and application to glioma.” Proc Natl Acad Sci, 104:20007-20012. [doi:10.1073/pnas.0710052104](https://doi.org/10.1073/pnas.0710052104)


Documentation
* [GISTIC homepage](https://software.broadinstitute.org/cancer/cga/gistic)


Important Notes
* Module Name: gistic (see [the modules page](/apps/modules.html) for more information)
* environment variables set
	+ GISTIC\_HOME* Example files in $GISTIC\_HOME/examplefiles



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session based on $GISTIC\_HOME/run\_gistic\_example:



```

[user@biowulf]$ **sinteractive**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load gistic**
[user@cn3144 ~]$ **cd /data/$USER**
[user@cn3144 user]$ **basedir=$PWD/example\_results**
[user@cn3144 user]$ **mkdir $basedir**
[user@cn3144 user]$ **segfile=$GISTIC\_HOME/examplefiles/segmentationfile.txt**
[user@cn3144 user]$ **markersfile=$GISTIC\_HOME/examplefiles/markersfile.txt**
[user@cn3144 user]$ **refgenefile=$GISTIC\_HOME/refgenefiles/hg16.mat**
[user@cn3144 user]$ **alf=$GISTIC\_HOME/examplefiles/arraylistfile.txt**
[user@cn3144 user]$ **cnvfile=$GISTIC\_HOME/examplefiles/cnvfile.txt**
[user@cn3144 user]$ **gistic2 \
 -b $basedir \
 -seg $segfile \
 -mk $markersfile \
 -refgene $refgenefile \
 -alf $alf \
 -cnv $cnvfile \
 -genegistic 1 \
 -smallmem 1 \
 -broad 1 \
 -brlen 0.5 \
 -conf 0.90 \
 -armpeel 1 \
 -savegene 1 \
 -gcm extreme**
[user@cn3144 user]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. gistic.sh). For example:



```

#!/bin/sh
set -e

module load gistic

basedir=$PWD/example_results
mkdir $basedir

segfile=$GISTIC_HOME/examplefiles/segmentationfile.txt
markersfile=$GISTIC_HOME/examplefiles/markersfile.txt
refgenefile=$GISTIC_HOME/refgenefiles/hg16.mat
alf=$GISTIC_HOME/examplefiles/arraylistfile.txt
cnvfile=$GISTIC_HOME/examplefiles/cnvfile.txt

## call script that sets MCR environment and calls GISTIC executable 
gistic2 \
 -b $basedir \
 -seg $segfile \
 -mk $markersfile \
 -refgene $refgenefile \
 -alf $alf \
 -cnv $cnvfile \
 -genegistic 1 \
 -smallmem 0 \
 -broad 1 \
 -brlen 0.5 \
 -conf 0.90 \
 -armpeel 1 \
 -savegene 1 \
 -gcm extreme

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] gistic.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. gistic.swarm). For example:



```

cd /data/$USER/gistic/set1 \
 && mkdir results \
 && gistic2 -b results -seg segfile.txt -mk markersfile.txt -refgene $GISTIC_HOME/refgenefiles/hg38.UCSC.add_miR.160920.refgene.mat -alf arraylistfile.txt -cnv cnvfile.txt -genegistic 1 -smallmem 0 -broad 1 -brlen 0.5 -conf 0.90 -armpeel 1 -savegene 1 -gcm extreme
cd /data/$USER/gistic/set2 \
 && mkdir results \
 && gistic2 -b results -seg segfile.txt -mk markersfile.txt -refgene $GISTIC_HOME/refgenefiles/hg38.UCSC.add_miR.160920.refgene.mat -alf arraylistfile.txt -cnv cnvfile.txt -genegistic 1 -smallmem 0 -broad 1 -brlen 0.5 -conf 0.90 -armpeel 1 -savegene 1 -gcm extreme
cd /data/$USER/gistic/set3 \
 && mkdir results \
 && gistic2 -b results -seg segfile.txt -mk markersfile.txt -refgene $GISTIC_HOME/refgenefiles/hg38.UCSC.add_miR.160920.refgene.mat -alf arraylistfile.txt -cnv cnvfile.txt -genegistic 1 -smallmem 0 -broad 1 -brlen 0.5 -conf 0.90 -armpeel 1 -savegene 1 -gcm extreme
cd /data/$USER/gistic/set4 \
 && mkdir results \
 && gistic2 -b results -seg segfile.txt -mk markersfile.txt -refgene $GISTIC_HOME/refgenefiles/hg38.UCSC.add_miR.160920.refgene.mat -alf arraylistfile.txt -cnv cnvfile.txt -genegistic 1 -smallmem 0 -broad 1 -brlen 0.5 -conf 0.90 -armpeel 1 -savegene 1 -gcm extreme

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f gistic.swarm [-g #] [-t #] --module gistic
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g #  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t #  Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module gistic  Loads the GISTIC module for each subjob in the swarm 
 | |
 | |
 | |








