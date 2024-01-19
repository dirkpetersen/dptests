

document.querySelector('title').textContent = 'ont-fast5-api on Biowulf';
ont-fast5-api on Biowulf


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



A set of tools to manipulate HDF5 files of the Oxford Nanopore .fast5 file format.
This module only provides access to the command line tools, not the python API.



Documentation
* Tools are documented on [GitHub](https://github.com/nanoporetech/ont_fast5_api)


Important Notes
* Module Name: ont-fast5-api (see [the modules page](/apps/modules.html) for more information)
* Tools are multthreaded. Please match the number of threads to the number of allocated CPUs.
* Example files in `$ONT_FAST5_API_TEST_DATA`



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --gres=lscratch:100 --cpus-per-task=2 --mem=10g**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144]$ **cd /lscratch/${SLURM\_JOB\_ID}**
[user@cn3144]$ **module load ont-fast5-api**
[user@cn3144]$ **cp -rL ${ONT\_FAST5\_API\_TEST\_DATA:-none}/\* .**
[user@cn3144]$ **ls**
Zymo-GridION-EVEN-BB-SN
[user@cn3144]$ **find Zymo-GridION-EVEN-BB-SN -name '\*.fast5' -printf '.' | wc -c**
160000
[user@cn3144]$ **du -sh Zymo-GridION-EVEN-BB-SN**
11G     Zymo-GridION-EVEN-BB-SN/
[user@cn3144]$ **single\_to\_multi\_fast5 -i Zymo-GridION-EVEN-BB-SN \
 -s multi\_fast5 -n 10000 -t $SLURM\_CPUS\_PER\_TASK --recursive -c vbz**
[user@cn3144]$ **ls -lh multi\_fast5**
total 6.6G
-rw-r--r-- 1 user group 408M May 26 10:42 batch_0.fast5
-rw-r--r-- 1 user group 425M May 26 10:49 batch_10.fast5
-rw-r--r-- 1 user group 425M May 26 10:49 batch_11.fast5
...
[user@cn3144]$ **du -sh multi\_fast5**
6.6G    multi_fast5

## copy results back to shared storage before exiting
[user@cn3144]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. ont\_fast5\_api.sh) like so:



```

#!/bin/bash
wd=$PWD
module load ont-fast5-api
cd /lscratch/${SLURM_JOB_ID} || exit 1
module load ont-fast5-api/3.3.0
cp -rL ${ONT_FAST5_API_TEST_DATA:-none}/* .
single_to_multi_fast5 -i Zymo-GridION-EVEN-BB-SN \
    -s multi_fast5 -n 10000 -t $SLURM_CPUS_PER_TASK --recursive -c vbz
mv multi_fast5 $wd

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=4 --mem=10g --gres=lscratch:25 ont_fast5_api.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. ont\_fast5\_api.swarm). For example:



```

multi_to_single_fast5 -i batch_0.fast5 -s batch_0 -t $SLURM_CPUS_PER_TASK
multi_to_single_fast5 -i batch_1.fast5 -s batch_1 -t $SLURM_CPUS_PER_TASK
multi_to_single_fast5 -i batch_2.fast5 -s batch_2 -t $SLURM_CPUS_PER_TASK

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f ont_fast5_api.swarm -g 10 -t 4 --module ont-fast5-api/3.3.0
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g #  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t #  Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module ont\_fast5\_api  Loads the ont\_fast5\_api module for each subjob in the swarm 
 | |
 | |
 | |








