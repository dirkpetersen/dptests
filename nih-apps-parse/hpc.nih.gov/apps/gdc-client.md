

document.querySelector('title').textContent = 'gdc-client on Biowulf';
gdc-client on Biowulf


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


Raw sequence data, stored as BAM files, makes up the bulk of data stored at the [NCI Genomic Data Commons (GDC)](https://gdc-portal.nci.nih.gov/). The size of a single file can vary greatly. Most BAM files stored in the GDC are in the 50 MB - 40 GB size range, with some of the whole genome BAM files reaching sizes of 200-300 GB.

The GDC Data Transfer Tool provides an optimized method of transferring data to and from the GDC, and enables resumption of interrupted transfers.


Documentation
* To see options for gdc-client, type at the prompt:
```
gdc-client --help
```
* GDC Data Transfer Tool: <https://docs.gdc.cancer.gov/Data_Transfer_Tool/Users_Guide/Getting_Started/>
* Genomics Data Commons: <https://gdc.cancer.gov/>
* Genomics Data Commons Data Portal: <https://portal.gdc.cancer.gov/>


Important Notes
* Module Name: gdc-client (see [the modules page](/apps/modules.html) for more information)
* Example files in $GDC\_EXAMPLES



The vast majority of data available from the GDC is access controlled. Users will need to first register and obtain an [authentication token](https://gdc-docs.nci.nih.gov/API/Users_Guide/Authentication_and_Authorization/) to access the controlled data.


Interactive job
[Interactive data transfers](/docs/userguide.html#int) are best run on Helix, the dedicated interactive data transfer system. 
Helix has a direct connection to the internet, and does not go through one of the HPC proxy servers. Sample session:



```


[user@helix ~]$ **module load gdc-client**

[user@helix ~]$ **gdc-client download 22a29915-6712-4f7a-8dba-985ae9a1f005**
100% [##################################################################################] Time: 0:00:02   1.86 MB/s
100% [##################################################################################] Time: 0:00:00 177.78 kB/s
Successfully downloaded: 1


```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. gdc-client.sh). For example:



```

#!/bin/bash
module load gdc-client
gdc-client download 22a29915-6712-4f7a-8dba-985ae9a1f005

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] gdc-client.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. gdc-client.swarm). For example:



```

gdc-client download cce8411d-dc96-5597-b198-fb47c5cf691c
gdc-client download 574b02d5-2de1-5aab-be8d-2c9d251dde9e
gdc-client download ad22c8a4-7767-5427-9271-5f4b506a124c
gdc-client download 4b8af859-b9cd-52b1-bc64-9bfa5d816a5d

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f gdc-client.swarm [-g #] [-t #] --module gdc-client --maxrunning 10
```

where


|  |  |  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module gdc-client Loads the gdc-client module for each subjob in the swarm 
 | --maxrunning 10 Only allow 10 simultaneous downloads at a time 
 | |
 | |
 | |
 | |


gdc-client pulls data from the GDC Data Portal, which can be overwhelmed by high numbers of simultaneous downloads, causing individual swarm subjobs to fail. It is best to include --maxrunning 10 to prevent this overload.








