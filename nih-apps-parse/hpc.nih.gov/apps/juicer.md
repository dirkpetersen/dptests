

document.querySelector('title').textContent = 'juicer on Biowulf';
juicer on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Batch job](#sbatch) 
[Interactive job](#int) 
[Swarm of jobs](#swarm) 
 |



Juicer is a system for analyzing loop-resolution Hi-C experiments. It was
developed in the Aiden lab at Baylor College of Medicine/Rice University. [ [Juicer website](http://aidenlab.org/juicer/docs.html)]


The juicer pipeline submits jobs to the cluster and then exits. It should
therefore be run on the login node. Individual tools can be run as batch jobs
or interactively.


Juicer depends on reference files (bwa index plus chromosome sizes file) and
restriction enzyme files that are part of the central install. You can build
additional references yourself or ask staff to generate them centrally. See
`$JUICER/references` and `$JUICER/restriction_sites` for
available reference data.


In addition to juicer, juicebox for visualizinig juicer contact maps
is also available as a separate module. Note, however, that juicebox
might be better run on your local desktop.


### References:


* Neva C. Durand, Muhammad S. Shamim, Ido Machol, Suhas S. P. Rao, 
 Miriam H. Huntley, Eric S. Lander, and Erez Lieberman Aiden. 
 *Juicer provides a one-click system for analyzing loop-resolution
 Hi-C experiments.*. Cell Systems 3, 2016.
 [PubMed](https://www.ncbi.nlm.nih.gov/pubmed/27467249)  | 
 [PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5846465/) |  
 [Journal](https://www.sciencedirect.com/science/article/pii/S2405471216302198?via%3Dihub)


 
Documentation
* Juicer [GitHub repo](https://github.com/theaidenlab/juicer)
* Juicer [wiki](https://github.com/theaidenlab/juicer/wiki)
* Juicer [Home](http://aidenlab.org/documentation.html)


Important Notes
* Module Names: juicer and juicebox (see 
 [the modules page](/apps/modules.html) for more information)
* Unusual environment variables set: `$JUICER`
* Reference data in `$JUICER/references` and 
 `$JUICER/restriction_sites`.
* Example data in `$JUICER_TEST_DATA`




Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
The juicer pipeline works by creating a series of batch jobs and submitting
them all at once using job dependencies. The main script is very light weight
and has to be run on the login node.


First create a directory for your Juicer run. A subdirectory called 'fastq'
within there should contain the fastq files. For example:



```

/data/user/juicer_test
+-- fastq
|   +-- reads_R1.fastq.gz
|   +-- reads_R2.fastq.gz

```

The fastq files must be called `filename_R1.fastq[.gz]` and
`filename_R2.fastq[.gz]`, as those names are built into the script. If
your fastq files have different names, you can rename them or create symlinks.



For your first run we recommend using the test data, which can be copied
from `$JUICER_TEST_DATA` and uncompressed:



```

[user@biowulf]$ **cd /data/user/juicer\_test**
[user@biowulf]$ **cp -r $JUICER\_TEST\_DATA/fastq .**
[user@biowulf]$ **ls -lh**
total 1.2G
-rw-r--r-- 1 user group 157M May  9 10:08 HIC003_S2_L001_R1_001.fastq.gz
-rw-r--r-- 1 user group 167M May  9 10:08 HIC003_S2_L001_R2_001.fastq.gz

```

 Juicer will create subdirectories aligned, HIC\_tmp,
debug, splits. The HIC\_tmp subdirectory will get
deleted at the end of the run.


By default, running `juicer.sh` with no options will use the hg19
reference file, and the DpnII restriction site map:



```

[user@biowulf]$ **module load juicer**
[user@biowulf]$ **juicer.sh**
Running juicer version 1.5.6                                                                                          
(-: Looking for fastq files...fastq files exist:                                 
-rw-r--r-- 1 user group 157M Feb 10  2017 /data/user/juicer_test/fastq/HIC003_S2_L001_R1_001.fastq.gz
-rw-r--r-- 1 user group 167M Feb 10  2017 /data/user/juicer_test/fastq/HIC003_S2_L001_R2_001.fastq.gz
(-: Aligning files matching /data/user/juicer_test/fastq/*_R*.fastq*
 in queue norm to genome hg19 with site file /usr/local/apps/juicer/juicer-1.5.6/restriction_sites/hg19_DpnII.txt
(-: Created /data/user/test_data/juicer/test_pipeline/splits and /data/user/test_data/juicer/test_pipeline/aligned
.                                                                                                                 
(-: Starting job to launch other jobs once splitting is complete
Submitted ligation counting job a1525867074HIC003_S2_L001_001.fastq_Count_Ligation
Submitted align1 job a1525867074_align1_HIC003_S2_L001_001.fastq.gz
Submitted align2 job a1525867074_align2_HIC003_S2_L001_001.fastq.gz
Submitted merge job a1525867074_merge_HIC003_S2_L001_001.fastq.gz
Submitted check job a1525867074_check
Submitted fragmerge job a1525867074_fragmerge
Submitted dedup_guard job in held state a1525867074_dedup_guard
Submitted dedup job a1525867074_dedup
Submitted post_dedup job a1525867074_post_dedup
Submitted dupcheck job a1525867074_dupcheck
Submitted stats job a1525867074_stats
Submitted hic job a1525867074_hic
Submitted hic job a1525867074_hic30
Submitted hiccups_wrap job a1525867074_hiccups_wrap
Submitted arrowhead job a1525867074_arrowhead_wrap
Submitted fincln job a1525867074_prep_done
(-: Finished adding all jobs... Now is a good time to get that cup of coffee..

```

 After the juicer.sh script exits, you should see a set of jobs in running
and queued state, with dependencies:



```

[user@biowulf]$ **squeue -u user**
 JOBID PARTITION     NAME   USER ST    TIME  NODES NODELIST(REASON)
466500      norm a1525867   user PD    0:00      1 (None)
466501      norm a1525867   user PD    0:00      1 (None)
466502      norm a1525867   user PD    0:00      1 (Dependency)
466504      norm a1525867   user PD    0:00      1 (Dependency)
466498      norm a1525867   user PD    0:00      1 (None)
466499      norm a1525867   user PD    0:00      1 (None)
466506      norm a1525867   user PD    0:00      1 (Dependency)
466507      norm a1525867   user PD    0:00      1 (Dependency)
466508      norm a1525867   user PD    0:00      1 (Dependency)
466509      norm a1525867   user PD    0:00      1 (Dependency)
466510      norm a1525867   user PD    0:00      1 (Dependency)
466511      norm a1525867   user PD    0:00      1 (Dependency)
466512      norm a1525867   user PD    0:00      1 (Dependency)
466513      norm a1525867   user PD    0:00      1 (Dependency)
466514      norm a1525867   user PD    0:00      1 (Dependency)
466503      norm a1525867   user PD    0:00      1 (Dependency)
466505      norm a1525867   user PD    0:00      1 (JobHeldUser)

```

You can follow the progress of the job by watching the jobs proceed, and by
examining the files in the debug subdirectory.


In addition to running the whole pipeline, individual batch jobs can also be run
as usual. For example, to use the `juicer_tools pre` command to
create a `.hic` format file from your own processed data create
a batch script similar to the following (juicertools.sh):



```

#! /bin/bash
module load juicer/1.6
juicer_tools pre $JUICER_TEST_DATA/input.txt.gz test.hic hg19

```

We have modified juicer\_tools to accept java options starting with -X. This
allows users to run with larger memory - for example `juicer_tools -Xmx48g ...`


Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=2 --mem=18g juicertools.sh
```

Some of the juicer\_tools require GPUs and therefore have to be submitted
to the gpu partition. For example (hiccups.sh):



```

#!/bin/bash

module load juicer
module load CUDA/8.0
juicer_tools  hiccups -m 500 -r 5000 -k KR -f 0.1 -p 4 \
    -i 10 -t 0.01,1.5,1.75,2 --ignore_sparsity test.hic output

# use juicer_tools48g if more memory is required

```

 Submit this job to a GPU node with:



```

[user@biowulf]$ sbatch  -p gpu --mem=18g  --gres=gpu:k80:1 hiccups.sh

```

Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
The juicer pipeline is not suitable for an interactive job. However,
individual tools can be run interactively. For this, allocate an [interactive session](/docs/userguide.html#int) and run the program.
Sample session:



```

[user@biowulf]$ **sinteractive --mem=18g --gres=gpu:k80:1**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144]$ **module load juicer**
[user@cn3144]$ **module load CUDA/8.0**
[user@cn3144]$ **juicer\_tools hiccups -m 500 -r 5000 -k KR -f 0.1 -p 4 \
 -i 10 -t 0.01,1.5,1.75,2 test.hic test\_output**

[user@cn3144]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf]$

```

In addition, an interactive session can also be used to run juicebox,
the visualization tool for HiC data. Note that X11 forwarding has to
be set up for this to work which will be the case if the sinteractive 
session is started from an NX session.



```

[user@biowulf]$ **sinteractive --mem=6g**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144]$ **module load juicebox**
[user@cn3144]$ **juicebox**

```

Should open the juicebox graphical user interface.


Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
The juicer pipeline isn't suitable for swarm. However, individual tools
can be run as a swarm. For this, create a swarmfile (e.g. juicer.swarm). For example:



```

juicer_tools pre sample1.txt.gz sample1.hic hg19
juicer_tools pre sample2.txt.gz sample2.hic hg19
juicer_tools pre sample3.txt.gz sample3.hic hg19

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f juicer.swarm -g 18 -t 1 --module juicer/1.5.6
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g #  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t #  Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module juicer  Loads the juicer module for each subjob in the swarm 
 | |
 | |
 | |








