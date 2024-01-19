

document.querySelector('title').textContent = 'Manta on Biowulf';
Manta on Biowulf


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


 Manta is a packaged used to discover structural variants and indels from
next generation sequencing data. It is optimized for rapid clinical
analysis, calling structural variants, medium-sized indels and large insertions.
Manta makes use of split read and paired end information and includes scoring
models optimized for germline analysis of diploid genomes and tumor-normal 
genome comparisons. Major use cases (as listed in the manta manual):



* Joint analysis of small sets of diploid individuals (where 'small' 
 means family-scale -- roughly 10 or fewer samples)
* Subtractive analysis of a matched tumor/normal sample pair
* Analysis of an individual tumor sample


There is also experimental RNA-Seq support.


### References:


* Chen et al. Manta: Rapid detection of structural variants and indels for clinical sequencing applications. 2015: http://dx.doi.org/10.1101/024232 [BioRxiv](http://biorxiv.org/content/early/2015/08/10/024232)


Documentation
* [Manta Main Site](https://github.com/Illumina/manta)


Important Notes
* Module Name: manta (see [the modules page](/apps/modules.html) for more information)
* Multithreaded/singlethreaded/MPI...
* environment variables set 
	+ $MANTA\_TEST\_DATA* Example files in $MANTA\_TEST\_DATA



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive -c 10 --mem 10g**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load manta**

[user@cn3144 ~]$ **configManta.py \**
  **--normalBam=${MANTA\_TEST\_DATA}/HCC1954.NORMAL.30x.compare.COST16011\_region.bam \**
  **--tumorBam=${MANTA\_TEST\_DATA}/G15512.HCC1954.1.COST16011\_region.bam \**
  **--referenceFasta=${MANTA\_TEST\_DATA}/Homo\_sapiens\_assembly19.COST16011\_region.fa \**
  **--region=8:107652000-107655000 \**
  **--region=11:94974000-94989000 \**
  **--candidateBins=4 --exome --runDir=./test**

[user@cn3144 ~]$ **tree test**
test
|-- [user   4.0K]  results
|   |-- [user   4.0K]  stats
|   `-- [user   4.0K]  variants
|-- [user   7.0K]  runWorkflow.py
|-- [user   3.0K]  runWorkflow.py.config.pickle
`-- [user   4.0K]  workspace

[user@cn3144 ~]$ **test/runWorkflow.py -m local -j 10 -g 10**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```

The workflow is executed by running the generated runWorkflow.py script. In our case, this is wrapped into a slurm batch script



Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. manta.sh). For example:



```

#!/bin/bash
module load manta || exit 1
test/runWorkflow.py -m local -j $SLURM_CPUS_PER_TASK -g $((SLURM_MEM_PER_NODE / 1024))

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=4 --mem=10g manta.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. manta.swarm). For example:



```

normal1_vs_tumor1/runWorkflow.py -m local -j $SLURM_CPUS_PER_TASK -g $((SLURM_MEM_PER_NODE / 1024))
normal2_vs_tumor2/runWorkflow.py -m local -j $SLURM_CPUS_PER_TASK -g $((SLURM_MEM_PER_NODE / 1024))
normal3_vs_tumor3/runWorkflow.py -m local -j $SLURM_CPUS_PER_TASK -g $((SLURM_MEM_PER_NODE / 1024))

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f manta.swarm -g 10 -t 10 --module manta
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module manta Loads the manta module for each subjob in the swarm 
 | |
 | |
 | |








