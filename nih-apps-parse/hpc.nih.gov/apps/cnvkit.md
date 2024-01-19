

document.querySelector('title').textContent = 'Cnvkit on HPC';
Cnvkit on HPC


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

  CNVkit is a Python library and command-line software toolkit to infer 
 and visualize copy number from targeted DNA sequencing data. It is designed 
 for use with hybrid capture, including both whole-exome and custom target 
 panels, and short-read sequencing platforms such as Illumina and Ion Torrent. 
 


Documentation * <http://cnvkit.readthedocs.io/en/stable/quickstart.html>



Important Notes
* Module Name: cnvkit (see [the modules 
 page](/apps/modules.html) for more information)
* Example files in /usr/local/apps/cnvkit/cnvkit-examples.

 To test cnvkit with the example files:
 
```

  $ cp -r /usr/local/apps/cnvkit/cnvkit-examples /data/$USER
  $ cd /data/$USER/cnvkit-examples
  $ sinteractive --mem=5g
  $ module load cnvkit
  $ make
  
```
* Reference data in 
 
```
/fdb/igenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa
/fdb/igenomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa
```
* Sequencing-accessible regions files are under

```
/usr/local/apps/cnvkit/data
	  
```





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

[user@cn3144 ~]$ **module load cnvkit
[user@cn3144 ~]$ cnvkit.py autobin input.bam**
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```




Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. cnvkit.sh). For example:



```

#!/bin/bash
set -e
module load cnvkit
cnvkit.py autobin input.bam
```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch cnvkit.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. cnvkit.swarm). For example:



```

cd dir1; cnvkit.py autobin input.bam
cd dir2; cnvkit.py autobin input.bam
cd dir3; cnvkit.py autobin input.bam

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f cnvkit.swarm [-g #] --module TEMPLATE
```

where
 

|  |  |
| --- | --- |
| -g *#*  | Number of Gigabytes of memory required for each process 
 (1 line in the swarm command file)  |
| --module  | Loads the module for each subjob in the swarm  |




