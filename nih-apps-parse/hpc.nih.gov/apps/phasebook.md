

document.querySelector('title').textContent = 'phasebook on Biowulf';
phasebook on Biowulf


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



phasebook is a novel approach for reconstructing the haplotypes of diploid genomes from long reads *de novo*, that is without the need for a reference genome.



### References:


* [Luo, Xiao, Xiongbin Kang, and Alexander Schönhuth. "phasebook: haplotype-aware de novo assembly of diploid genomes from long reads." *Genome biology* 22.1 (2021): 1-26.](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02512-x)


Documentation
* [phasebook on GitHub](https://github.com/phasebook/phasebook)


Important Notes
* Module Name: phasebook (see [the modules page](/apps/modules.html) for more information)
 * Use the $SLURM\_CPUS\_PER\_TASK variable to set the appropriate number of threads. (See example below.)
 * Environment variables set 
	+ $PHASEBOOK\_TESTDATA points to a small test data set* The online documentation suggest calling the main script like so: python phasebook.py. On our system it is incorrect to use the python prefix with this script due to the way it is installed.



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf ~]$ **sinteractive -c8 --mem=4g --gres=lscratch:10**
salloc: Pending job allocation 33141417
salloc: job 33141417 queued and waiting for resources
salloc: job 33141417 has been allocated resources
salloc: Granted job allocation 33141417
salloc: Waiting for resource configuration
salloc: Nodes cn0881 are ready for job
srun: error: x11: no local DISPLAY defined, skipping
error: unable to open file /tmp/slurm-spank-x11.33141417.0
slurmstepd: error: x11: unable to read DISPLAY value

[user@cn0881 ~]$ **module load phasebook**
[+] Loading phasebook  1.0.0  on cn0881
[+] Loading singularity  3.8.5-1  on cn0881

[user@cn0881 ~]$ **cd /lscratch/$SLURM\_JOB\_ID**

[user@cn0881 33141417]$ **cp -r $PHASEBOOK\_TESTDATA .**

[user@cn0881 33141417]$ **cd TESTDATA**

[user@cn0881 TESTDATA]$ **phasebook.py -i reads.fa -t $SLURM\_CPUS\_PER\_TASK -p hifi -g small -x**
use preset parameters...
2022-02-25 13:20:59,226 - /opt/phasebook/scripts/phasebook.py[line:262] - INFO: splitting input fastx file into 1 subfiles...
2022-02-25 13:20:59,361 - /opt/phasebook/scripts/phasebook.py[line:270] - INFO: splitting finished.
[...snip]
start polishing...
2022-02-25 13:22:26,529 - /opt/phasebook/scripts/phasebook.py[line:374] - INFO: All has been finished successfully.

2022-02-25 13:22:26,529 - /opt/phasebook/scripts/phasebook.py[line:375] - INFO: The final output haplotype aware contigs are here: ./contigs.fa

2022-02-25 13:22:26,529 - /opt/phasebook/scripts/phasebook.py[line:376] - INFO: Thank you for using phasebook!

[user@cn0881 TESTDATA]$ **ls**
1.split_fastx  3.cluster        5.polish          clustered_reads.list  phasebook.log  reference.fa
2.overlap      4.asm_supereads  all.supereads.fa  contigs.fa            reads.fa

[user@cn0881 TESTDATA]$ **exit**
exit
salloc: Relinquishing job allocation 33141417

[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. phasebook.sh). For example:



```

#!/bin/bash
set -e
mkdir -p /data/${USER}/phasebook/${SLURM_JOB_ID}
cd /data/${USER}/phasebook/${SLURM_JOB_ID}
module load phasebook
phasebook.py -i /path/to/reads.fa -t $SLURM_CPUS_PER_TASK -p hifi -g small -x

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] phasebook.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. phasebook.swarm). For example:



```

mkdir -p /data/${USER}/phasebook/${SLURM_JOB_ID}; cd !$; phasebook.py -i /path/to/reads1.fa -t $SLURM_CPUS_PER_TASK -p hifi -g small -x
mkdir -p /data/${USER}/phasebook/${SLURM_JOB_ID}; cd !$; phasebook.py -i /path/to/reads2.fa -t $SLURM_CPUS_PER_TASK -p hifi -g small -x
mkdir -p /data/${USER}/phasebook/${SLURM_JOB_ID}; cd !$; phasebook.py -i /path/to/reads3.fa -t $SLURM_CPUS_PER_TASK -p hifi -g small -x
mkdir -p /data/${USER}/phasebook/${SLURM_JOB_ID}; cd !$; phasebook.py -i /path/to/reads4.fa -t $SLURM_CPUS_PER_TASK -p hifi -g small -x

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f phasebook.swarm [-g #] [-t #] --module phasebook
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module phasebook Loads the phasebook module for each subjob in the swarm 
 | |
 | |
 | |








