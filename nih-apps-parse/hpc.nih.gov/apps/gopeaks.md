

document.querySelector('title').textContent = "gopeaks";
gopeaks on Biowulf


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



From the gopeaks documentation:




> 
>  GoPeaks is a peak caller designed for CUT&TAG/CUT&RUN sequencing
>  data. GoPeaks by default works best with narrow peaks such as H3K4me3 and
>  transcription factors. However, broad epigenetic marks like H3K27Ac/H3K4me1
>  require different the step, slide, and minwidth parameters.
> 


### References:


* W. M. Yashar  *et al*
[**GoPeaks: histone modification peak calling for CUT&Tag**](https://pubmed.ncbi.nlm.nih.gov/35788238/) 
 *Genome Biol. 2022 Jul 4;23(1):14*


Documentation
* [Gopeaks on GitHub](https://github.com/maxsonBraunLab/gopeaks)


Important Notes
* Module Name: gopeaks (see [the modules page](/apps/modules.html) for more information)
* gopeaks can take advantage of more than one CPU though it does not appear to be
 configurable via well known environment variables or command line arguments. Based
 on testing an allocation of 4-6 CPUs is appropriate
* Example files in $GOPEAKS\_TEST\_DATA* Reference data in /fdb/gopeaks/



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**) for running the peakfinder with a control and generating a summary plot with deeptools:



```

[user@biowulf]$ **sinteractive --cpus-per-task=6 --gres=lscratch:20 --mem=30g**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load gopeaks**

[user@cn3144 ~]$ **cd /lscratch/$SLURM\_JOB\_ID**
[user@cn3144 ~]$ **cp $GOPEAKS\_TEST\_DATA/\* .**
[user@cn3144 ~]$ **ls -lh**

[user@cn3144 ~]$ **gopeaks -b GSE190793\_Kasumi\_cutrun.bam -c GSE190793\_Kasumi\_IgG.bam --mdist 1000 --prefix Kasumi\_cnr**
Reading chromsizes from bam header...
nTests: 431593
nzSignals: 4.5114732e+07
nzBins: 8298786
n: 7.566515e+06
p: 7.184688090883924e-07
mu: 5.962418894299423
var: 5.962414610487421

[user@cn3144 ~]$ **cat Kasumi\_cnr\_gopeaks.json**
{
        "gopeaks_version": "1.0.0",
        "date": "2023-05-12 12:40:59 PM",
        "elapsed": "7m19.045741504s",
        "prefix": "Kasumi_cnr",
        "command": "gopeaks -b GSE190793_Kasumi_cutrun.bam -c GSE190793_Kasumi_IgG.bam --mdist 1000 --verbose --prefix=Kasumi_cnr",
        "peak_counts": 29393
}

[user@cn3144 ~]$ # create a summary graph for peaks centered on the middle of the peak interval +- 1000
[user@cn3144 ~]$ # i.e. not gene annotation

[user@cn3144 ~]$ **module load deeptools**
[user@cn3144 ~]$ **bamCoverage -p6 -b GSE190793\_Kasumi\_cutrun.bam -o GSE190793\_Kasumi\_cutrun.bw**
[user@cn3144 ~]$ **bamCoverage -p6 -b GSE190793\_Kasumi\_IgG.bam -o GSE190793\_Kasumi\_IgG.bw**
[user@cn3144 ~]$ **computeMatrix reference-point -R Kasumi\_cnr\_peaks.bed -a 1000 -b 1000 --referencePoint center \
 -S GSE190793\_Kasumi\_cutrun.bw GSE190793\_Kasumi\_IgG.bw \
 --sortRegions descend --samplesLabel 'Cut&Run' 'IgG' -p6 -o cutrun\_matrix**
[user@cn3144 ~]$ **plotHeatmap -m cutrun\_matrix -o cutrun.png --averageTypeSummaryPlot mean --colorMap GnBu**
[user@cn3144 ~]$ **cp Kasumi\_cnr\_\* \*.bw cutrun.png /data/$USER/my\_working\_directory**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


![](/images/gopeaks_heatmap.png)


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. gopeaks.sh). For example:



```

#!/bin/bash
set -e
module load gopeaks/1.0.0
cp $GOPEAKS_TEST_DATA/* .
gopeaks -b GSE190793_Kasumi_cutrun.bam -c GSE190793_Kasumi_IgG.bam --mdist 1000 --prefix Kasumi_cnr

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=6 --mem=30g gopeaks.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. gopeaks.swarm). For example:



```

gopeaks -b replicate1.bam -c control.bam --mdist 1000 --prefix replicate1_gopeaks
gopeaks -b replicate2.bam -c control.bam --mdist 1000 --prefix replicate2_gopeaks

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f gopeaks.swarm -g 30 -t 6 --module gopeaks
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module gopeaks Loads the gopeaks module for each subjob in the swarm 
 | |
 | |
 | |








