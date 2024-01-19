

document.querySelector('title').textContent = 'bonito on Biowulf';
bonito on Biowulf


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




 A PyTorch Basecaller for Oxford Nanopore Reads. According to ONT this is a research
 release
 



> 
>  provided as technology demonstrators to provide early access to
>  features or stimulate Community development of tools. Support for this
>  software will be minimal and is only provided directly by the
>  developers.
>  



Documentation
* bonito on [GitHub](https://github.com/nanoporetech/bonito)


Important Notes
* Module Name: bonito (see [the modules page](/apps/modules.html) for more information)
* This software is GPU accelerated
* Example files in `$BONITO_TEST_DATA`


### Models available


Use



```

$ bonito download --models --show

```

to see all available models. Note that models are included in the install
already. *fast*: fast model *hac*: high accuracy *sup*:
super high accuracy


Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Note that 0.3.6
runs on p100 GPUs but >=0.5.0 requires v100 or newer. Sample session:



```

[user@biowulf]$ **sinteractive --gres=lscratch:50,gpu:v100x:1 --mem=12g --cpus-per-task=6**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144]$ **cd /lscratch/$SLURM\_JOB\_ID**
[user@cn3144]$ **module load bonito/0.7.2**
[user@cn3144]$ **cp -rL ${BONITO\_TEST\_DATA:-none}/\* .**
[user@cn3144]$ **ls -lh**
total 4.0K
drwxr-xr-x 3 user group 4.0K Feb  8 11:07 Zymo-GridION-EVEN-BB-SN
[user@cn3144]$ **find Zymo-GridION-EVEN-BB-SN -name '\*.fast5' -printf '.' | wc -c**
160000
[user@cn3144]$ ### basecalling command for bonito 0.3.6
[user@cn3144]$ **bonito basecaller --fastq --recursive \
 --device cuda dna\_r9.4.1 Zymo-GridION-EVEN-BB-SN > reads.fastq**
> loading model
> calling: 20829 reads [49:08,  7.08 reads/s]
...
[user@cn3144]$ ### basecalling command for bonito >=0.5.0
[user@cn3144]$ **bonito basecaller --recursive --device cuda \
 dna\_r9.4.1\_e8\_hac@v3.3 Zymo-GridION-EVEN-BB-SN | gzip -c - > reads.fastq.gz**
> loading model dna_r9.4.1_e8.1_hac@v3.3
> completed reads: 160000
> duration: 0:24:34
> samples per second 4.7E+06
> done

[user@cn3144]$ ### basecalling command for bonito >=0.7.1
[user@cn3144]$ **bonito basecaller --recursive --device cuda \
 dna\_r9.4.1\_e8\_hac@v3.3 Zymo-GridION-EVEN-BB-SN | gzip -c - > reads.fastq.gz**
> reading fast5
> outputting unaligned fastq
> loading model dna_r9.4.1_e8_hac@v3.3
> completed reads: 160000
> duration: 0:23:39
> samples per second 4.9E+06
> done

[user@cn3144]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf]$


```

Current versions of bonito achieved ~110 reads/s on a V100X GPU with the hac
model. The basecaller does not scale to more than 1 GPU



Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. bonito.sh), which uses the input file 'bonito.in'. For example:



```

#!/bin/bash
wd=$PWD

module load bonito/0.7.2 || exit 1
cd /lscratch/$SLURM_JOB_ID || exit 1
cp -rL ${BONITO_TEST_DATA:-none}/* .
bonito basecaller --recursive --device cuda \
    dna_r9.4.1_e8_hac@v3.3 Zymo-GridION-EVEN-BB-SN | gzip -c - > reads.fastq.gz
mv reads.fastq $wd

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=6 --mem=12g --gres=gpu:v100x:1,lscratch:50 bonito.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. bonito.swarm). For example:



```

bonito basecaller --fastq --recursive --device cuda dna_r9.4.1_e8_hac@v3.3 run1 > reads1.fastq
bonito basecaller --fastq --recursive --device cuda dna_r9.4.1_e8_hac@v3.3 run2 > reads2.fastq

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f bonito.swarm -g 12 -t 6 --gres=gpu:v100x:1 --module bonito/0.7.2
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g #  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t #  Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module bonito  Loads the bonito module for each subjob in the swarm
 | |
 | |
 | |






