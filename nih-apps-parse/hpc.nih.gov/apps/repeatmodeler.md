

document.querySelector('title').textContent = 'RepeatModeler on Biowulf';
RepeatModeler on Biowulf


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



RepeatModeler is a de novo transposable element (TE) family identification and modeling package. At the heart of RepeatModeler are three de-novo repeat finding programs ( RECON, RepeatScout and LtrHarvest/Ltr\_retriever ) which employ complementary computational methods for identifying repeat element boundaries and family relationships from sequence data.



Documentation
* [RepeatModeler Main Site](http://www.repeatmasker.org/RepeatModeler/)


Important Notes
* Module Name: repeatmodeler (see [the modules page](/apps/modules.html) for more information)
 * It's important to set the -pa # option/argument pair properly. From the RepeatModeler help section:
   

Specify the number of parallel search jobs to run. RMBlast jobs will
 use 4 cores each and ABBlast jobs will use a single core each. i.e.
 on a machine with 12 cores and running with RMBlast you would use
 -pa 3 to fully utilize the machine.



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf ~]$ **sinteractive -c4 --mem=4g --gres=lscratch:10**
salloc.exe: Pending job allocation 47680219
salloc.exe: job 47680219 queued and waiting for resources
salloc.exe: job 47680219 has been allocated resources
salloc.exe: Granted job allocation 47680219
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn0857 are ready for job
srun: error: x11: no local DISPLAY defined, skipping

[user@cn0857 ~]$ **cd /lscratch/$SLURM\_JOB\_ID**

[user@cn0857 47680219]$ **cp /fdb/app\_testdata/fasta/R64-1-1.cdna\_nc.fa .**

[user@cn0857 47680219]$ **module load repeatmodeler**
[+] Loading repeatmodeler  2.0.1  on cn0857
[+] Loading singularity  3.5.2  on cn0857

[user@cn0857 47680219]$ **BuildDatabase -name test R64-1-1.cdna\_nc.fa**
Building database test:
  Reading R64-1-1.cdna_nc.fa...
Number of sequences (bp) added to database: 7126 ( 9153986 bp )

[user@cn0857 47680219]$ **RepeatModeler -database test -pa 1 -LTRStruct >& run.out #runs for ~20m**

[user@cn0857 47680219]$ **tail run.out**
The results have been saved to:
  test-families.fa  - Consensus sequences for each family identified.
  test-families.stk - Seed alignments for each family identified.

The RepeatModeler stockholm file is formatted so that it can
easily be submitted to the Dfam database.  Please consider contributing
curated families to this open database and be a part of this growing
community resource.  For more information contact help@dfam.org.

[user@cn0857 47680219]$ **exit**
exit
salloc.exe: Relinquishing job allocation 47680219

[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. repeatmodeler.sh). For example:



```

#!/bin/bash
set -e
module load repeatmodeler
RepeatModeler -database test -pa 1 -LTRStruct >& run.out

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] repeatmodeler.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. repeatmodeler.swarm). For example:



```

RepeatModeler -database a -pa 1 -LTRStruct >& run_a.out
RepeatModeler -database b -pa 1 -LTRStruct >& run_b.out
RepeatModeler -database c -pa 1 -LTRStruct >& run_c.out
RepeatModeler -database d -pa 1 -LTRStruct >& run_d.out

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f repeatmodeler.swarm [-g #] [-t #] --module repeatmodeler
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module repeatmodeler Loads the repeatmodeler module for each subjob in the swarm 
 | |
 | |
 | |








