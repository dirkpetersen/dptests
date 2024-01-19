

document.querySelector('title').textContent = 'Bamreadcount on HPC';
Bamreadcount on HPC


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


 [Bam-readcount](https://github.com/genome/bam-readcount) generates 
 metrics at single nucleotide positions.  

 There are number of metrics generated which can be useful for filtering 
 out false positive calls.



Documentation
* <https://github.com/genome/bam-readcount>



Important Notes
* Module Name: bamreadcount (see [the 
 modules page](/apps/modules.html) for more information)
* cram-v0.0.1 version support CRAM. Example files in $BAMREADCOUNT\_TEST\_DATA.





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

[user@cn3144 ~]$ **module load bamreadcount**
[user@cn3144 ~]$ **cp -rp /usr/local/apps/bamreadcount/0.8.0/test-data /data/$USER**[user@cn3144 ~]$ **cd /data/$USER/test-data**[user@cn3144 ~]$ **bam-readcount -f ref.fa test.bam**
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
module load bamreadcount
bam-readcount -f ref.fa test.bam
```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch batch.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. job.swarm). For example:



```

cd dir1; bam-readcount -f ref.fa test.bam
cd dir2; bam-readcount -f ref.fa test.bam
cd dir3; bam-readcount -f ref.fa test.bam

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f job.swarm --module bamreadcount
```

where
 

|  |  |
| --- | --- |
| -g *#*  | Number of Gigabytes of memory required for each process 
 (1 line in the swarm command file)  |
| --module  | Loads the module for each subjob in the swarm  |




