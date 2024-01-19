

document.querySelector('title').textContent = 'Subread on Biowulf';
Subread on Biowulf


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


Subread package: high-performance read alignment, quantification and mutation discovery
 The Subread package comprises a suite of software programs for processing next-gen sequencing read data including:


* Subread: an accurate and efficient aligner for 
 mapping both genomic DNA-seq reads and RNA-seq reads (for the purpose 
 of expression analysis).
* Subjunc: an RNA-seq aligner suitable for all purposes 
 of RNA-seq analyses.
* featureCounts: a highly efficient and accurate 
 read summarization program.
* exactSNP: a SNP caller that discovers SNPs by testing signals against local background noises.


Subread can be run multi-threaded using -T flag on biowulf. See example below.
### References:

 * <https://www.ncbi.nlm.nih.gov/pubmed/23558742>
* <http://www.ncbi.nlm.nih.gov/pubmed/24227677>


Documentation * <http://subread.sourceforge.net/>


Important Notes * Module Name: subread (see [the 
 modules page](/apps/modules.html) for more information) 
 * Subread can be run multi-threaded using -T flag on biowulf. See example 
 below.



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load subread**
[user@cn3144 ~]$ **cd /data/$USER/dir**
[user@cn3144 ~]$ **subread-align -i indexfile -r inputfile -o outputfile**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit). Create a batch input file (e.g. subread.sh). For example:



```

#!/bin/bash
set -e
module load subread
cd /data/$USER/dir
subread-align -i indexfile -r inputfile -o outputfile
```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
$ sbatch subread.sh
```

To run multi-threaded subread on biowulf:



```

#!/bin/bash

module load subread
cd /data/$USER/
subread-align -T $SLURM_CPUS_PER_TASK -i indexfile -r inputfile -o output
```


Submit the script:
 
```
$ sbatch --cpus-per-task=4 jobscript
```


 --cpus-per-task: allocate 4 cpus. This number will be assigned to $SLURM\_CPUS\_PER\_TASK automatically 
 For more memory requirement (default 4gb), use --mem flag: 
 
```
$ sbatch --cpus-per-task=4 --mem=10g jobscript
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. subread.swarm). For example:



```

cd dir1;subread-align -i indexfile -r inputfile -o outputfile
cd dir2;subread-align -i indexfile -r inputfile -o outputfile
cd dir3;subread-align -i indexfile -r inputfile -o outputfile
cd dir4;subread-align -i indexfile -r inputfile -o outputfile

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f subread.swarm [-g #] [-t #] --module subread
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module subread Loads the subread module for each subjob in the swarm  | |
 | |
 | |




















