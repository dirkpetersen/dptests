

document.querySelector('title').textContent = "pod5";
pod5 on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
 |



From the documentation




> 
>  POD5 is a file format for storing nanopore dna data in an easily accessible way. The format is able to be written in a streaming manner which allows a sequencing instrument to directly write the format.
>  Data in POD5 is stored using Apache Arrow, allowing users to consume data in many languages using standard tools.
> 



Documentation
* pod5 repo on [GitHub](https://github.com/nanoporetech/pod5-file-format)
* pod5 [Manual](https://pod5-file-format.readthedocs.io/en/latest/)


Important Notes
* Module Name: pod5 (see [the modules page](/apps/modules.html) for more information)
* This tool is multithreaded. Please match threads to allocated CPUs
* Environment variables set 
	+ $POD5\_TEST\_DATA* Example files in $POD5\_TEST\_DATA



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive --cpus-per-task=6 --mem=12g --gres=lscratch:150**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load pod5**

[user@cn3144 ~]$ **cd /lscratch/$SLURM\_JOB\_ID**
[user@cn3144 ~]$ **cp -rL ${POD5\_TEST\_DATA:-none} input**
[user@cn3144 ~]$ **pod5 convert fast5 --threads $SLURM\_CPUS\_PER\_TASK --output output.pod5 --recursive input**
[user@cn3144 ~]$ **du -sh input**
46G     input
[user@cn3144 ~]$ **ls -lh output.pod5**
-rw-r--r-- 1 user group  38G Jun  2 14:29 output.pod5
[user@cn3144 ~]$ **pod5 view --output summary.tsv output.pod5**
[user@cn3144 ~]$ **head summary.tsv**
[user@cn3144 ~]$ **pod5 inspect read output.pod5 0001297c-4c07-438e-a29b-6da3b0ad1260**
read_id: 0001297c-4c07-438e-a29b-6da3b0ad1260
read_number:    11392
start_sample:   180540114
median_before:  220.86135864257812
channel data:
        channel: 284
        well: 1
        pore_type: not_set
end reason:
        name: unknown
        forced: False
...

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```

Note that `pod5 convert fast5` requires multi-fast5 input files



Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. pod5.sh). For example:



```

#!/bin/bash
set -e
module load pod5

cd /lscratch/$SLURM_JOB_ID
mkdir output
cp -rL ${POD5_TEST_DATA:-none} input
pod5 convert fast5 --threads $SLURM_CPUS_PER_TASK --output output input/* 
cp 

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=4 --mem=10g --gres=lscratch:150 pod5.sh
```







