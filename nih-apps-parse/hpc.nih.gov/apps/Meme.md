

document.querySelector('title').textContent = 'Meme on Biowulf';
Meme on Biowulf


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



The MEME Suite allows you to:
* discover motifs using MEME, DREME (DNA only) or GLAM2 on groups of related DNA or protein sequences,
* search sequence databases with motifs using MAST, FIMO, MCAST or GLAM2SCAN,
* compare a motif to all motifs in a database of motifs,
* associate motifs with Gene Ontology terms via their putative target genes, and
* analyse motif enrichment using SpaMo or CentriMo.


The Meme Suite was developed at U. Queensland and U. Washington. 
[Meme website](http://meme-suite.org).
Meme is cpu-intensive for large numbers of sequences or long sequences and 
scales well to 128 cores. 

Meme motif and GoMo databases are available in /fdb/meme/

[![meme](/images/meme-programs.jpg)](http://meme-suite.org) 




Documentation
* Type 'meme' or any of the Meme Suite programs with no parameters on the command line to see a list
of all available options and more information.
* [Meme Suite documentation](http://meme-suite.org)


Important Notes
* Module Name: meme (see [the modules page](/apps/modules.html) for more information)
* Uses MPI for parallelization. 
* Your input database should consist of a file containing sequences in fasta format. There are several example files in /usr/local/apps/meme/examples* Reference data in /fdb/meme/* **Maxsize parameter:** The maximum dataset size in characters. Determine the number of characters in your dataset by typing 'wc -c filename'. e.g.

```

[user@biowulf mydir]$ wc -c mini-drosoph.s 
506016 mini-drosoph.s

```


For this dataset, the maxsize parameter has to be set to greater than 506,016, so we will use 600000 as the Maxsize parameter. See example under 'Batch job' below.



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --ntasks=4 --ntasks-per-core=1**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$  **cd /data/$USER**

[user@cn3144 ~]$ **module load meme**
+] Loading meme  4.12.0  on cn3144
[+] Loading openmpi 2.1.1  for GCC 4.8.5

[user@cn3144]$ **cp /usr/local/apps/meme/examples/protease-seqs /data/$USER**

[user@cn3144]$ **meme -text protease-seqs -p $SLURM\_NTASKS > protease.meme.out**
IInitializing the motif probability tables for 2 to 7 sites...
nsites = 7
Done initializing.
SEEDS: highwater mark: seq 6 pos 300
BALANCE: samples 7 residues 1750 nodes 1 residues/node 1750

seqs=     7, min= 185, max=  300, total=     1750

motif=1
SEED WIDTHS: 8 11 15 21 29 41 50
em: w=  50, psites=   7, iter=   0

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).

In the example below, we are using the file 'mini-drosoph.s' as the input. This file can be copied from /usr/local/apps/meme/examples. The maxsize parameter will be set to 
600000 as described in the [Important Notes](#notes") section above.

Create a batch input file (e.g. meme.sh). For example:



```

#!/bin/bash
set -e

cd /data/$USER/
module load meme
meme mini-drosoph.s  -oc meme_out -maxsize 600000 -p $SLURM_NTASKS

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --ntasks=28  --ntasks-per-core=1 --exclusive meme.sh
```


where
* The meme command in the batch script does not include 'mpirun', or specify the nodes to be used. Meme uses OpenMPI that has built-in Slurm support, so it will automatically get the node/cpu list from Slurm. 
* The number of MPI processes will be pulled from the variable $SLURM\_NTASKS. (-p $SLURM\_NTASKS in the meme command in the batch script).
* --ntasks=28 : run 28 MPI tasks 
* --ntasks-per-core=1 : run only one MPI task per physical core (i.e. don't use hyperthreading)
* --exclusive : allocate the node exclusively



This job will run on the norm (default) partition, which is limited to single-node jobs. It will utilize all the cores on a 28-core norm node. Meme scales well 

and large meme jobs (maxsize ~500,000) can be submitted on up to 512 cores. A multinode job can be submitted with:

```
sbatch --partition=multinode --ntasks=512 --ntasks-per-core=1  meme.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources. 

You would submit a swarm of Meme jobs if you have several input files.

Meme is an MPI program which uses OpenMPI libraries. OpenMPI on Biowulf is
built with Slurm support. An MPI program runs a specified number of MPI processes or 'tasks'. The user specifies the number of tasks with '--ntasks=#' on the sbatch command line, and the
OpenMPI program automatically gets this number from Slurm and starts up the appropriate number of tasks.

Swarm is intended for single-threaded and multi-threaded applications.
When you use the '-t #' (threads per process) flag to swarm, it sets up
subjobs with $SLURM\_CPUS\_PER\_TASK=# and allocates # cpus on a single
node for each subjob. The Meme MPI program sees this as a single 'task'
with #threads, and not as # tasks, and will complain that there are not enough slots available
for the MPI processes.

Thus, it is important to add the flag --sbatch '--ntasks=# when submitting a swarm of Meme jobs. You should also use '--ntasks-per-core=1' as most MPI applications
run with greater efficiency with only one MPI task on each physical core. 

Create a swarmfile (e.g. Meme.swarm). For example:



```

meme query1.fa -oc query1.out -maxsize 10000000 -p $SLURM_NTASKS
meme query2.fa -oc query2out -maxsize 10000000 -p $SLURM_NTASKS
meme query3.fa -oc query3.out -maxsize 10000000 -p $SLURM_NTASKS

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f swarm.cmd -g 20 --sbatch '--ntasks=4 --ntasks-per-core=1' --module=meme
```





























