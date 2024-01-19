

document.querySelector('title').textContent = "hifiasm";
hifiasm on Biowulf


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



 Hifiasm is a fast haplotype-resolved de novo assembler initially designed for PacBio HiFi reads.
 Its latest release supports telomere-to-telomere assembly by utilizing ultralong Oxford Nanopore reads.
 It can produce better haplotype-resolved assemblies when given parental short reads or Hi-C data.



### References:


* Cheng, H., Concepcion, G.T., Feng, X., Zhang, H., Li H. (2021)
 [**Haplotype-resolved de novo assembly using phased assembly graphs with hifiasm.**](https://doi.org/10.1038/s41592-020-01056-5)
 Nat Methods, 18:170-175. doi:10.1038/s41592-020-01056-5
* Cheng, H., Jarvis, E.D., Fedrigo, O., Koepfli, K.P., Urban, L., Gemmell, N.J., Li, H. (2022)
 [**Haplotype-resolved assembly of diploid genomes without parental data.**](https://doi.org/10.1038/s41587-022-01261-x)
 Nature Biotechnology, 40:1332–1335. doi:10.1038/s41587-022-01261-x


Documentation
* [hifiasm Main Site](https://hifiasm.readthedocs.io/en/latest/index.html)
* [hifiasm on GitHub](https://github.com/chhylp123/hifiasm)


Important Notes
* Module Name: hifiasm (see [the modules page](/apps/modules.html) for more information)
* Multithreaded. Set number of threads with the -t flag



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive --cpus-per-task 32**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load hifiasm**

[user@cn3144 ~]$ **hifiasm -t $SLURM\_CPUS\_PER\_TASK -o NA24385\_180901\_0113427.asm /fdb/app\_testdata/fastq/Homo\_sapiens/HG002\_NA24385\_son-PacBio\_CCS\_15kb/m54238\_180901\_011437.Q20.fastq.gz**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. hifiasm.sh). For example:



```

#!/bin/bash
set -e
module load hifiasm
hifiasm -t $SLURM_CPUS_PER_TASK -o NA24385_180901_0113427.asm /fdb/app_testdata/fastq/Homo_sapiens/HG002_NA24385_son-PacBio_CCS_15kb/m54238_180901_011437.Q20.fastq.gz

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=# [--mem=#] hifiasm.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. hifiasm.swarm). For example:



```

hifiasm -t $SLURM_CPUS_PER_TASK -o sample1.asm sample1.fq.gz
hifiasm -t $SLURM_CPUS_PER_TASK -o sample2.asm sample2.fq.gz
hifiasm -t $SLURM_CPUS_PER_TASK -o sample3.asm sample3.fq.gz
hifiasm -t $SLURM_CPUS_PER_TASK -o sample4.asm sample4.fq.gz

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f hifiasm.swarm [-g #] -t # --module hifiasm
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module hifiasm Loads the hifiasm module for each subjob in the swarm 
 | |
 | |
 | |








