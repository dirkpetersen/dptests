

document.querySelector('title').textContent = 'RSeQC on Biowulf';
RSeQC on Biowulf


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


 RSeQC package provides a number of useful modules that can comprehensively evaluate high throughput sequence data especially RNA-seq data. Some basic modules quickly inspect sequence quality, nucleotide composition bias, PCR bias and GC bias, while RNA-seq specific modules evaluate sequencing saturation, mapped reads distribution, coverage uniformity, strand specificity, transcript level RNA integrity etc.

### References:

 * <https://academic.oup.com/bioinformatics/article/28/16/2184/325191>
* <https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-0922-z>


Documentation * <http://rseqc.sourceforge.net/>


Important Notes * Module Name: rseqc (see [the modules 
 page](/apps/modules.html) for more information)



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

[user@cn3144 ~]$ **module load rseqc**[user@cn3144 ~]$ **bam\_stat.py -i input.bam $> output**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit). Create a batch input file (e.g. rseqc.sh). For example:



```

#!/bin/bash
set -e
module load rseqc
bam_stat.py -i input.bam > output
```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.

 
```
sbatch [--mem=#] rseqc.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources. Create a swarmfile (e.g. rseqc.swarm). For example:



```
cd dir1;bam_stat.py -i input.bam > output
cd dir2;bam_stat.py -i input.bam > output
cd dir3;bam_stat.py -i input.bam > output
```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f rseqc.swarm [-g #] --module rseqc
```

where
 

|  |  |  |  |
| --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in 
 the swarm command file) 
 | --module rseqc Loads the TEMPLATE module for each subjob in the swarm 
  | |
 | |










