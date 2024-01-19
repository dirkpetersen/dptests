

document.querySelector('title').textContent = 'Misopy on HPC';
Misopy on HPC


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
[Swarm of jobs](#swarm) 
[Parallel jobs](#parallel) 
 |


 MISO (Mixture-of-Isoforms) is a probabilistic framework that quantitates 
 the expression level of alternatively spliced genes from RNA-Seq data, and 
 identifies differentially regulated isoforms or exons across samples. By 
 modeling the generative process by which reads are produced from isoforms 
 in RNA-Seq, the MISO model uses Bayesian inference to compute the probability 
 that a read originated from a particular isoform.


### References:


* The MISO framework is described in Katz et. al., [Analysis 
 and design of RNA sequencing experiments for identifying isoform regulation](http://www.nature.com/nmeth/journal/v7/n12/full/nmeth.1528.html). 
 *Nature Methods* (2010).



Documentation
* <http://miso.readthedocs.org/en/fastmiso/>



Important Notes
* Module Name: misopy (see [the modules 
 page](/apps/modules.html) for more information)
* Multithreaded/parallel
* miso\_setting.txt
 
  

 The default miso setting files can be copied from */usr/local/apps/misopy/miso\_settings.txt* 
 and modified. It looks like this: 

```
[data]
filter_results = True
min_event_reads = 20

[cluster]
cluster_command = sbatch --mem=20g --cpus-per-task=8 --time=99:00:00

[sampler]
burn_in = 500
lag = 10
num_iters = 5000
num_chains = 6
num_processors = 8

```

 Miso can be run in two modes: multi-threaded mode and parallel mode. Each 
 mode use/ignore different directives in the miso setting file. 
 


Multi-threaded mode will be used when running miso without '--use-cluster' 
 flag.   

 The line 'num\_processors' will be used and 'cluster\_command" will 
 be ignored. The default threads are 8. 


Parallel mode will be used when running miso with both '--use-cluster' 
 and '--chunk-jobs' flags on biowulf. In this mode, depending on 
 the number assigned to '--chunk-jobs=#', many jobs will be submitted 
 to the cluster and each job will use one thread no matter what is 
 assigned to 'num\_processors' line since this line will be ignored. 
 Events will beatch and each batch will be submitted to a job. The 
 smaller the # is, the more number of jobs will be created. The default 
 memory is 20gb and 99 hours of walltime for each job. If more memory 
 or walltime is needed, copy/modify the setting file, include '--settings-filename=The\_Full\_Path\_to\_Settingfile' 
 flag in miso command, then submit miso job using sbatch. 


If '--settings-filename=.....' flag is not specified in miso command, the 
 default miso setting file will be used.






 


Interactive job
[Interactive jobs](/docs/userguide.html#int) 
 should be used for debugging, graphics, or applications that cannot be run 
 as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) 
 and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=20g --cpus-per-task=8**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load misopy**
[user@cn3144 ~]$ **miso --run ./indexed ./accepted\_hits.bam --output-dir miso\_out --read-len 76**
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```

 


Batch job
Most jobs should be run as [batch 
 jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. batch.sh). For example:



```

#!/bin/bash
set -e
module load misopy
cd dir
miso --run ./indexed ./accepted_hits.bam --output-dir miso_out --read-len 76
```

Submit this job using the Slurm [sbatch](/docs/userguide.html) 
 command.



```
sbatch --cpus-per-task=8 --mem=20g batch.sh
```


Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is 
 an easy way to submit a set of independent commands requiring identical 
 resources.
Create a swarmfile (e.g. job.swarm). For example:



```

cd dir1; miso --run ./indexed ./accepted_hits.bam --output-dir miso_out --read-len 76
cd dir2; miso --run ./indexed ./accepted_hits.bam --output-dir miso_out --read-len 76
cd dir3; miso --run ./indexed ./accepted_hits.bam --output-dir miso_out --read-len 76

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f job.swarm -g 20 -t 8 --module misopy
```

 where 
 

|  |  |
| --- | --- |
| -g *#*  | Number of Gigabytes of memory required for each process 
 (1 line in the swarm command file)  |
| -t *#* | Number of threads/CPUs required for each process (1 line in the swarm 
 command file).  |
| --module  | Loads the module for each subjob in the swarm  |




Submit a Parallel Job
Refer to [Notes section](#Notes) regarding miso setting file. 1. Create a script file similar to the lines below.
 



```
#!/bin/bash
module load misopy
cd /data/$USER/dir
miso --run \
./indexed \
./accepted_hits.bam \
--output-dir miso_out_cluster \
--read-len 76 \
--use-cluster \
--chunk-jobs=1000
```

2. Submit the script on biowulf:



```
$ sbatch jobscript
```

In this example, events will be split into multiple chunks, each chunk 
 contains 1000 events, and each chunk will be submitted as a single job to 
 the cluster. By default each job will run on 1 core and 20gb or memory. 
 The default walltime is 99 hours.


For more memory or walltime, copy and modify the miso\_setting.txt file and 
 add '--settings-filename=/data/$USER/..../miso\_settings\_cluster.txt' to the 
 miso command. 






