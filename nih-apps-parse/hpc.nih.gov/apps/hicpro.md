

document.querySelector('title').textContent = 'HiC-Pro on HPC';
HiC-Pro on HPC


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

  HiC-Pro was designed to process Hi-C data, from raw fastq files (paired-end 
 Illumina data) to the normalized contact maps. Since version 2.7.0, HiC-Pro 
 supports the main Hi-C protocols, including digestion protocols as well 
 as protocols that do not require restriction enzyme such as DNase Hi-C. 
 In practice, HiC-Pro can be used to process dilution Hi-C, in situ Hi-C, 
 DNase Hi-C, Micro-C, capture-C, capture Hi-C or HiChip data.  

 The pipeline is flexible, scalable and optimized. It can operate either 
 on a single laptop or on a computational cluster. HiC-Pro is sequential 
 and each step of the workflow can be run independantly.  

 HiC-Pro includes a fast implementatation of the iterative correction method 
 (see the iced python library for more information).


### References:

 * <http://www.genomebiology.com/2015/16/1/259>


Documentation
* <http://nservant.github.io/HiC-Pro/>



Important Notes
* Module Name: hicpro (see [the modules 
 page](/apps/modules.html) for more information)
* Test files: /usr/local/apps/hicpro/test
* The --cpur-per-task should match the N\_CPU set in the config.txt. For example, if it's 'N\_CPU = 4' in the config.txt, then one should use --cpus-per-task=4 when submit the job.





Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --cpus-per-task=4 --mem=40g**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load hicpro**
[user@cn3144 ~]$ **cp -r /usr/local/apps/hicpro/test /data/$USER/hicpro/test**
[user@cn3144 ~]$ **cd /data/$USER/hicpro/test**
[user@cn3144 ~]$ **tar xvfz HiCPro\_testdata.tar.gz**

# Modify the config.txt and change 'username@hpc.nih.gov' to your email address

[user@cn3144 ~]$ **HiC-Pro -i test\_data -o output -c config.txt**

# Running the aligment step of HiC-Pro in parallel.
[user@cn3144 ~]$ **module load hicpro/3.1.0\_v2**

# First, split reads in small chunks, by default --nreads would be 2,000,000:
[user@cn3144 ~]$ **split\_reads.py \
 --results\_folder split\_test/dixon\_2M\_split \
 --nreads 50000 \
 ./dixon\_2M/SRR400264\_00\_R1.fastq.gz** 
[user@cn3144 ~]$ **split\_reads.py \
 --results\_folder split\_test/dixon\_2M\_split \
 --nreads 50000 \
 ./dixon\_2M/SRR400264\_00\_R2.fastq.gz** 

# Then, generate two sbatch jobs:
[user@cn3144 ~]$ **HiC-Pro \
 -i split\_test \
 -o split\_out \
 -c config.txt -p

Run HiC-Pro 3.1.0 parallel mode
The following command will launch the parallel workflow through 5 torque jobs:
sbatch HiCPro\_step1\_apptest1.sh
The following command will merge the processed data and run the remaining steps per sample:
sbatch HiCPro\_step2\_apptest1.sh** 

# Last, submit sbatch jobs:
[user@cn3144 ~]$ **cd split\_out**
[user@cn3144 ~]$ **sbatch --dependency=afterok:$(sbatch HiCPro\_step1\_apptest1.sh) HiCPro\_step2\_apptest1.sh** 
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```




Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. sbatch.sh). For example:



```

#!/bin/bash
set -e
module load hicpro
cd /data/$USER/hicpro/test
HiC-Pro -i test_data -o output -c config.txt
```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=4 --mem=40g sbatch.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. job.swarm). For example:



```

cd dir1; HiC-Pro -i test_data -o output -c config.txt
cd dir2; HiC-Pro -i test_data -o output -c config.txt
cd dir3; HiC-Pro -i test_data -o output -c config.txt
```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f job.swarm -g 40 -t 4 --module hicpro
```

where
 

|  |  |
| --- | --- |
| -g *#*  | Number of Gigabytes of memory required for each process (1 line in the swarm command file)
  |
| -t *#* | Number of threads/CPUs required for each process (1 line in the swarm command file).
  |
| --module  | Loads the module for each subjob in the swarm  |




