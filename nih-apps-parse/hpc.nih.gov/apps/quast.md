

document.querySelector('title').textContent = 'QUAST on Biowulf';
QUAST on Biowulf


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



 QUAST stands for QUality ASsessment Tool. The tool evaluates genome assemblies by computing various metrics. The package includes the general QUAST tool for genome assemblies, MetaQUAST, the extension for metagenomic datasets, and Icarus, interactive visualizer for these tools.



### References:


* Alla Mikheenko, Andrey Prjibelski, Vladislav Saveliev, Dmitry Antipov, Alexey Gurevich,
 Versatile genome assembly evaluation with QUAST-LG,
 Bioinformatics (2018) 34 (13): i142-i150. doi: [10.1093/bioinformatics/bty266](https://doi.org/10.1093/bioinformatics/bty266)
 First published online: June 27, 2018
* Alla Mikheenko, Gleb Valin, Andrey Prjibelski, Vladislav Saveliev, Alexey Gurevich,
 Icarus: visualizer for de novo assembly evaluation,
 Bioinformatics (2016). doi: [10.1093/bioinformatics/btw379](https://doi.org/10.1093/bioinformatics/btw379)
 First published online: July 4, 2016
* Alla Mikheenko, Vladislav Saveliev, Alexey Gurevich,
 MetaQUAST: evaluation of metagenome assemblies,
 Bioinformatics (2016) 32 (7): 1088-1090. doi: [10.1093/bioinformatics/btv697](https://doi.org/10.1093/bioinformatics/btv697)
 First published online: November 26, 2015
* Alexey Gurevich, Vladislav Saveliev, Nikolay Vyahhi and Glenn Tesler,
 QUAST: quality assessment tool for genome assemblies,
 Bioinformatics (2013) 29 (8): 1072-1075. doi: [10.1093/bioinformatics/btt086](https://doi.org/10.1093/bioinformatics/btt086)
 First published online: February 19, 2013


Documentation
* [QUAST Main Site](http://quast.sourceforge.net)
* [QUAST manual](http://quast.bioinf.spbau.ru/manual.html)


Important Notes
* Module Name: quast (see [the modules page](/apps/modules.html) for more information)
 * This application produces HTML reports. You can use [hpcdrive](/docs/hpcdrive.html) to view these reports on your local workstation.
* Environment variables set 
	+ QUAST\_HOME* Example files in $QUAST\_HOME/test\_data



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive --gres lscratch:10**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load quast**
[+] Loading quast, version 5.0.2...
[user@cn3144 ~]$ **cd /lscratch/$SLURM\_JOB\_ID**
[user@cn3144 ~]$ **quast.py --threads 2 \
 $QUAST\_HOME/test\_data/contigs\_{1,2}.fasta \
 -R $QUAST\_HOME/test\_data/reference.fasta.gz \
 -G $QUAST\_HOME/test\_data/genes.gff** 

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. quast.sh). For example:



```

#!/bin/bash
set -e
module load quast

test -n "$SLURM_CPUS_PER_TASK" || SLURM_CPUS_PER_TASK=2

quast.py --threads $SLURM_CPUS_PER_TASK \
    $QUAST_HOME/test_data/contigs_{1,2}.fasta \
    -R $QUAST_HOME/test_data/reference.fasta.gz \
    -G $QUAST_HOME/test_data/genes.gff

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] quast.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. quast.swarm). For example:



```

quast.py $QUAST_HOME/test_data/contigs_{1,2}.fasta \
         -R $QUAST_HOME/test_data/reference.fasta.gz \
         -G $QUAST_HOME/test_data/genes.gff
quast.py $QUAST_HOME/test_data/contigs_{1,2}.fasta \
         -R $QUAST_HOME/test_data/reference.fasta.gz \
         -G $QUAST_HOME/test_data/genes.gff
quast.py $QUAST_HOME/test_data/contigs_{1,2}.fasta \
         -R $QUAST_HOME/test_data/reference.fasta.gz \
         -G $QUAST_HOME/test_data/genes.gff

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f quast.swarm [-g #] -t # --module quast
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module quast Loads the QUAST module for each subjob in the swarm 
 | |
 | |
 | |








