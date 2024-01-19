

document.querySelector('title').textContent = 'Bsmap on HPC';
Bsmap on HPC


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

  BSMAP is a short reads mapping software for bisulfite sequencing reads. Bisulfite 
 treatment converts unmethylated Cytosines into Uracils (sequenced as Thymine) 
 and leave methylated Cytosines unchanged, hence provides a way to study 
 DNA cytosine methylation at single nucleotide resolution. BSMAP aligns the 
 Ts in the reads to both Cs and Ts in the reference. 


RRBSMAP is a specifically designed version of BSMAP for reduced representation 
 bisulfite sequencing (RRBS), it indexes the genome only on the enzyme digestion 
 sites and therefore guarantees all reads were mapped to digestion sites, 
 and greatly reduces the CPU/memory usage. Since BSMAP-2.0, RRBSMAP has been 
 merged into BSMAP.



### References:

 * <https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-10-232>


Documentation * <https://github.com/genome-vendor/bsmap>



Important Notes * Module Name: bsmap (see [the modules 
 page](/apps/modules.html) for more information)
* Multithreaded





Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --cpus-per-task=4**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load bsmap**
[user@cn3144 ~]$ **bsmap -a infile -d ref.fa -o out.bam -p $SLURM\_CPUS\_PER\_TASK**

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
module load bsmap
bsmap -a infile -d ref.fa -o out.bam -p $SLURM_CPUS_PER_TASK
```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=4 batch.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. job.swarm). For example:



```

cd dir1; bsmap -a infile -d ref.fa -o out.bam -p $SLURM_CPUS_PER_TASK
cd dir2; bsmap -a infile -d ref.fa -o out.bam -p $SLURM_CPUS_PER_TASK
cd dir3; bsmap -a infile -d ref.fa -o out.bam -p $SLURM_CPUS_PER_TASK

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f job.swarm -t 4 --module bsmap
```

where
 

|  |  |
| --- | --- |
| -t *#* | Number of threads/CPUs required for each process (1 line 
 in the swarm command file).  |
| --module  | Loads the module for each subjob in the swarm  |




