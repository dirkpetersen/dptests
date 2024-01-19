

document.querySelector('title').textContent = 'Pindel on HPC';
Pindel on HPC


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

  Pindel can detect breakpoints of large deletions, medium sized insertions, 
 inversions, tandem duplications and other structural variants at single-based 
 resolution from next-gen sequence data. It uses a pattern growth approach 
 to identify the breakpoints of these variants from paired-end short reads.


### References:

 * <http://bioinformatics.oxfordjournals.org/content/25/21/2865.full.pdf>


Documentation * <http://gmt.genome.wustl.edu/packages/pindel/>



Important Notes * Module Name: pindel (see [the modules 
 page](/apps/modules.html) for more information)
* environment variables set to /data/$USER/pindeltempdir
	+ TMPDIR
* Example files in /usr/local/apps/pindel/TEST\_DATA


 


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

[user@cn3144 ~]$ **module load pindel**
[user@cn3144 ~]$ **cd /data/$USER/; cp ${PINDEL\_TEST\_DATA:-none}/\* .**
[user@cn3144 ~]$ **pindel -i simulated\_config.txt -f simulated\_reference.fa -o outfile -c ALL**

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
module load pindel
pindel -i simulated_config.txt -f simulated_reference.fa -o outfile -c ALL
```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--mem=#] batch.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. job.swarm). For example:



```
cd dir1; pindel -i simulated_config.txt -f simulated_reference.fa -o outfile -c ALL
cd dir2; pindel -i simulated_config.txt -f simulated_reference.fa -o outfile -c ALL

cd dir3; pindel -i simulated_config.txt -f simulated_reference.fa -o outfile -c ALL
```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f job.swarm [-g #] --module pindel
```

where
 

|  |  |
| --- | --- |
| -g *#*  | Number of Gigabytes of memory required for each process 
 (1 line in the swarm command file)  |
| --module  | Loads the module for each subjob in the swarm  |




