

document.querySelector('title').textContent = 'Csvkit on HPC';
Csvkit on HPC


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

  Csvkit is a suite of command-line tools for converting to and working with 
 CSV, the king of tabular file formats.

 ### 


Documentation * <https://csvkit.readthedocs.io/en/1.0.3/>



Important Notes * Module Name: csvkit (see [the modules 
 page](/apps/modules.html) for more information)
* example : /usr/local/apps/csvkit/csvkit\_tutorial/ne\_1033.data.xlsx





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

[user@cn3144 ~]$ **module load csvkit**
[user@cn3144 ~]$ **cd /data/$USER/dir**
[user@cn3144 ~]$ **curl -L -O https://github.com/onyxfish/csvkit/raw/master/examples/realdata/ne\_1033\_data.xlsx**
[user@cn3144 ~]$ **in2csv /data/$USER/dir/ne\_1033\_data.xlsx**

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
module load csvkit
in2csv exampleInfile
```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--mem=#] batch.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. job.swarm). For example:



```

cd dir1; in2csv infile
cd dir2; in2csv infile
cd dir3; in2csv infile

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f job.swarm [-g #] --module csvkit
```

where
 

|  |  |
| --- | --- |
| -g *#*  | Number of Gigabytes of memory required for each process 
 (1 line in the swarm command file)  |
| --module  | Loads the module for each subjob in the swarm  |




