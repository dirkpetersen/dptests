

document.querySelector('title').textContent = 'bcl-convert on Biowulf';
bcl-convert on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
 |



Illumina's bcl-convert is the (future) successor to bcl2fastq. The application converts Binary Base Call (BCL)
files produced by Illumina sequencing systems to FASTQ files. bcl-convert also provides adapter handling
(through masking and trimming) and UMI trimming and produces metric outputs.




The current setup of bcl-convert on Biowulf requires **at least 32 CPUs**, and preferrably a whole node,
to run without overloading a compute node. Additionally certain options have been preset and setting
them will cause an error. See [Important Notes](#notes) below!



Documentation
* [bcl-convert Documentation](https://support-docs.illumina.com/SW/BCL_Convert/Content/SW/FrontPages/BCL_Convert.htm)
* [Upgrading from bcl2fastq to BCL Convert](https://support.illumina.com/bulletins/2020/10/upgrading-from-bcl2fastq-to-bcl-convert.html)


Important Notes
* Module Name: bcl-convert (see [the modules page](/apps/modules.html) for more information)
* Multithreaded. **Requires a minimum of 32 CPUs**.
* bcl-convert will write logs to /data/$USER/.bclconvert\_logs or to lscratch, if available


### bcl-convert Options


Do **NOT** set the following options when running bcl-convert:
 * --bcl-num-parallel-tiles
* --bcl-num-conversion-threads
* --bcl-num-compression-threads
* --bcl-num-decompression-threads





### sbatch/sinteractive Options


You **MUST** set the following sbatch/sinteractive options as described below.




|  |  |
| --- | --- |
| **Option** | **Explanation/Howto** |
| --exclusive | The node must be allocated exclusively, else your bcl-convert process will overload CPUs and be inefficient/run slower. |
| --constraint | The number of CPUs on the allocated node must be known, so that bcl-convert will run the correct number of threads. To determine
 this, use the freen command to find the different types of nodes and select one type. (example below) |
| --cpus-per-task |  Must be set to the number of CPUs on the node type you are requesting. |
| --gres=lscratch | *Optional*, bcl-convert will write temporary logs in lscratch. Additionally using lscratch to write output may be beneficial. See example session below |
| --mem | *Optional*, set to all the available memory on the type of node you are requesting. |




Example session to choose parameters

```

biowulf% **freen**
                                                    .......Per-Node Resources......
Partition    FreeNds      FreeCPUs           Cores  CPUs  GPUs    Mem   Disk Features
-------------------------------------------------------------------------------------------------------
norm*         0 / 695    7926 / 38920         28    56         247g   400g cpu56,core28,g256,ssd400,x2695,ibfdr
norm*         0 / 521    2806 / 29168         28    56         247g   800g cpu56,core28,g256,ssd800,x2680,ibfdr
norm*         0 / 7       144 / 392           28    56         247g  2400g cpu56,core28,g256,ssd2400,x2680,ibfdr
norm*       397 / 519   13290 / 16608         16    32         121g   800g cpu32,core16,g128,ssd800,x2650,10g
[...]

```

freen reports that there are 'cpu32' (32 CPUs), or 'cpu56' (56 CPUs) nodes available.   

Thus, to submit to a 32-cpu node (121 GB of RAM), your sbatch or sinteractive command would have the parameters: 

```
--exclusive --constraint=cpu32 --cpus-per-task=32 --mem=121g --gres=lscratch:400
```

Alternatively, if you want a 56-CPU node, you would use:

```
--exclusive --constraint=cpu56 --cpus-per-task=56 --mem=247g --gres=lscratch:400
```

Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive --exclusive --constraint=cpu32 --cpus-per-task=32 --mem=64G --gres=lscratch:400**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **cd /lscratch/$SLURM\_JOBID**

[user@cn3144 ~]$ **mkdir sample\_bclconvert\_output**

[user@cn3144 ~]$ **module load bcl-convert**

[user@cn3144 ~]$ **bcl-convert --bcl-input-directory /data/$USER/sample-run \
 --output-directory sample\_bclconvert\_output**
Index Read 2 is marked as Reverse Complement in RunInfo.xml: The barcode and UMI outputs will be output in Reverse Complement of Sample Sheet inputs.
Sample sheet being processed by common lib? Yes
SampleSheet Settings:
  AdapterRead1 = CAAGCAGAAGACGGCATACGAGAT
  AdapterRead2 = CAAGCAGAAGACGGCATACGAGAT
  FastqCompressionFormat = gzip
  SoftwareVersion = 3.7.4

shared-thread-linux-native-asio output is disabled
bcl-convert Version 00.000.000.3.9.3
Copyright (c) 2014-2018 Illumina, Inc.
...
[user@cn3144 ~]$ **mv sample\_bclconvert\_output /data/$USER/**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. bcl-convert.sh). For example:



```

#!/bin/bash
set -e
mkdir -p /lscratch/$SLURM_JOBID/sample-output
module load bcl-convert
bcl-convert --bcl-input-directory sample-run --output-directory /lscratch/$SLURM_JOBID/sample-output
mv /lscratch/$SLURM_JOBID/sample-output /data/$USER/

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --exclusive --constraint=cpu32 --cpus-per-task=32 --mem=64G --gres=lscratch:400 bcl-convert.sh
```









