

document.querySelector('title').textContent = 'Ceas on Biowulf';
Ceas on Biowulf


|  |
| --- |
| 
Quick Links
[References](#ref)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
[Swarm of jobs](#swarm) 
 |



CEAS is a tool designed to characterize genome-wide protein-DNA interaction patterns from ChIP-chip and ChIP-Seq of both sharp and broad binding factors. It provides statistics on ChIP enrichment at important genome features such as specific chromosome, promoters, gene bodies, or exons, and infers genes most likely to be regulated by a binding factor. CEAS also enables biologists to visualize the average ChIP enrichment signals over specific genomic features, allowing continuous and broad ChIP enrichment to be perceived which might be too subtle to detect from ChIP peaks alone. 



References
* [Shin, Hyunjin, et al. "CEAS: cis-regulatory element annotation system." *Bioinformatics* 25.19 (2009): 2605-2606.](https://academic.oup.com/bioinformatics/article/25/19/2605/182052?login=true)


Important Notes
* Module Name: ceas (see [the modules page](/apps/modules.html) for more information)
* Singlethreaded app
* Reference data in /fdb/CEAS/



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

[user@cn3144 ~]$**module load ceas**

[user@cn3144 ~]$ **cd /data/$USER/ceas**

[user@cn3144 ~]$ **cp /fdb/CEAS/\* .**

[user@cn3144 ~]$ **ceas --name=H3K36me3\_ceas --pf-res=20 --gn-group-names='Top 10%,Bottom 10%' \**
  **-g /fdb/CEAS/hg18.refGene -b H3K36me3\_MACS\_pval1e-5\_peaks.bed -w H3K36me3.wig**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. ceas.sh). For example:



```

#!/bin/bash
module load ceas
cd /data/$USER/ceas
ceas --name=H3K36me3_ceas --pf-res=20 --gn-group-names='Top 10%,Bottom 10%'\
  -g /fdb/CEAS/hg18.refGene -b H3K36me3_MACS_pval1e-5_peaks.bed -w H3K36me3.wig

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --mem=#10g ceas.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. ceas.swarm). For example:



```

  cd /data/$USER/dir1; ceas commands
  cd /data/$USER/dir2; ceas commands
  cd /data/$USER/dir3; ceas commands

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f ceas.swarm -g 10 --module ceas
```

where


|  |  |  |  |
| --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | --module ceas Loads the ceas module for each subjob in the swarm 
 | |
 | |








