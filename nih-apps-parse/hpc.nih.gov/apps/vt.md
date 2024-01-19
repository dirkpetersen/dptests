

document.querySelector('title').textContent = 'VT on HPC';
VT on HPC


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

 [VT](http://genome.sph.umich.edu/wiki/Vt#Installation) is a variant 
 tool set that discovers short variants from Next Generation Sequencing data.


### References:

 * <http://bioinformatics.oxfordjournals.org/content/31/13/2202>


Documentation * <https://genome.sph.umich.edu/wiki/Vt#How_to_cite_vt.3F>



Important Notes * Module Name: vt (see [the modules 
 page](/apps/modules.html) for more information)
* Example files in /usr/local/apps/vt/version/test





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

[user@cn3144 ~]$ **module load vt**
[user@cn3144 ~]$ **vt normalize IN.vcf -r ref.fa -o out.vcf**

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
module load vt
vt normalize IN.vcf -r ref.fa -o out.vcf
```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch batch.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. job.swarm). For example:



```

cd dir1; vt normalize IN.vcf -r ref.fa -o out.vcf
cd dri2; vt normalize IN.vcf -r ref.fa -o out.vcf
cd dir3; vt normalize IN.vcf -r ref.fa -o out.vcf

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f job.swarm [-g #] [-t #] --module vt
```

where
 

|  |  |
| --- | --- |
| -g *#*  | Number of Gigabytes of memory required for each process (1 line in the swarm command file)
  |
| -t *#* | Number of threads/CPUs required for each process (1 line in the swarm command file).
  |
| --module  | Loads the module for each subjob in the swarm  |




