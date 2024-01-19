

document.querySelector('title').textContent = 'RGT on HPC';
RGT on HPC


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


 HINT (Hmm-based IdeNtification of Transcription factor footprints) integrates 
 both DNase I hypersensitivity and histone modifications for the detection 
 of open chromatin regions and active binding sites. Within transcription 
 factor binding sites, there is a specific grammar of DNase I digestion and 
 histone marks. The authors have therefore devised a multivariate HMM to 
 model this regulatory grammar by simultaneous analysis of DNase-seq and 
 the ChIP-seq profiles of histone modifications on a genome-wide level. The 
 HMM has as input a normalized and a slope signal of DNase-seq and one of 
 the histone marks. It can therefore detect the increase, top and decrease 
 regions of either histone modification and DNase signals. The genomic regions 
 annotated with the HMM state are considered predictions and represent likely 
 binding sites within that cell context. For benchmarking data of main publication 
 please visit authors lab's [website](http://costalab.org/publications-2/hint-bc/ "HINT on CostaLab"). 
 


### References:


* <http://www.nature.com/nmeth/journal/vaop/ncurrent/full/nmeth.3772.html>


Documentation
* <http://www.regulatory-genomics.org/rgt/rgt-data-folder/>



Important Notes
* Module Name: rgt (see [the modules 
 page](/apps/modules.html) for more information)
* Example files in /usr/local/apps/rgt/DendriticCells.tar.gz and /usr/local/apps/rgt/HINT\_DNaseTest.tar.gz
* RGT data files in /fdb/rgt/rgtdata. The variable RGTDATA is set automatically when the rgt module is loaded





Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=10g**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load rgt**[user@cn3144 ~]$ **rgt-hint footprinting --dnase-seq DNase.bam DNasePeaks.bed**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```




Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. batch.sh). For example:



```

#!/bin/bash
#SBATCH --job-name="script"
#SBATCH --mail-type=BEGIN,END

cp -p /usr/local/apps/rgt/DendriticCells.tar.gz /data/$USER/rgt
tar xvfz /data/$USER/rgt/DendriticCells
cd /data/$USER/rgt/DendriticCells 

module load rgt/0.13.2

rgt-hint footprinting --atac-seq --paired-end --organism=mm10 --output-location=./ --output-prefix=pDC pDC.bam pDC_peaks.narrowPeak

rgt-hint footprinting --atac-seq --paired-end --organism=mm10 --output-location=./ --output-prefix=cDC1 cDC1.bam cDC1_peaks.narrowPeak 

rgt-hint tracks --bc --bigWig --organism=mm10 cDC1.bam cDC1_peaks.narrowPeak  --output-prefix=cDC1_BC
rgt-hint tracks --bc --bigWig --organism=mm10 pDC.bam pDC_peaks.narrowPeak  --output-prefix=pDC_BC

rgt-motifanalysis matching --organism=mm10 --input-files pDC.bed cDC1.bed
rgt-hint differential --organism=mm10 --bc --nc $SLURM_CPUS_PER_TASK --mpbs-files=./match/cDC1_mpbs.bed,./match/pDC_mpbs.bed --reads-files=cDC1.bam,pDC.bam --conditions=cDC1,pDC --output-location=cDC1_pDC

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command. The $SLURM\_CPUS\_PER\_TASK in the script will be replaced with the actual number of cpus automatically.



```
sbatch --cpus-per-task=16 --mem=25g --time=10:00:00 script
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. job.swarm). For example:



```

cd dir1; rgt-hint footprinting --dnase-seq DNase.bam DNasePeaks.bed
cd dir2; rgt-hint footprinting --dnase-seq DNase.bam DNasePeaks.bed
cd dir3; rgt-hint footprinting --dnase-seq DNase.bam DNasePeaks.bed

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f job.swarm [-g #] [-t #] --module rgt
```

where
 

|  |  |
| --- | --- |
| -g *#*  | Number of Gigabytes of memory required for each process 
 (1 line in the swarm command file)  |
| -t *#*  | Number of cpus required for each process
 (1 line in the swarm command file)  |
| --module  | Loads the module for each subjob in the swarm  |




