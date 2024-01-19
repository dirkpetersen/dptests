

document.querySelector('title').textContent = 'Conpair on HPC';
Conpair on HPC


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


Conpair: concordance and contamination estimator for tumor–normal pairs 
 


Conpair is a fast and robust method dedicated for human tumor-normal studies 
 to perform concordance verification (= samples coming from the same individual), 
 as well as cross-individual contamination level estimation in whole-genome 
 and whole-exome sequencing experiments. Importantly, our method of estimating 
 contamination in the tumor samples is not affected by copy number changes 
 and is able to detect contamination levels as low as 0.1%.



Documentation
* <https://github.com/nygenome/Conpair>



Important Notes
* Module Name: conpair (see [the 
 modules page](/apps/modules.html) for more information)
* Example files in /usr/local/apps/conpair/version/data/example





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

[user@cn3144 ~]$ **module load conpair**
[user@cn3144 ~]$ **cp -rp /usr/local/apps/conpair/10102016/data/example /data/$USER**[user@cn3144 ~]$ **cd /data/$USER/example/pileup**[user@cn3144 ~]$ **verify\_concordance.py -N NA12878\_normal40x.gatk.pileup.10lines.txt -T NA12878\_tumor80x.gatk.pileup.txt**
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```




Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. batch.sh). For example:



```

#!/bin/bash
set -e
module load conpair
verify_concordance.py -N NA12878_normal40x.gatk.pileup.10lines.txt -T NA12878_tumor80x.gatk.pileup.txt
```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch batch.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. job.swarm). For example:



```

cd dir1; verify_concordance.py -N NA12878_normal40x.gatk.pileup.10lines.txt -T NA12878_tumor80x.gatk.pileup.txt
cd dir2; verify_concordance.py -N NA12878_normal40x.gatk.pileup.10lines.txt -T NA12878_tumor80x.gatk.pileup.txt
cd dir3; verify_concordance.py -N NA12878_normal40x.gatk.pileup.10lines.txt -T NA12878_tumor80x.gatk.pileup.txt

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f job.swarm --module conpair
```

where
 

|  |  |
| --- | --- |
| -g *#*  | Number of Gigabytes of memory required for each process 
 (1 line in the swarm command file)  |
| --module  | Loads the module for each subjob in the swarm  |




