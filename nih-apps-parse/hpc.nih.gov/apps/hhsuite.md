

document.querySelector('title').textContent = 'hhsuite on Biowulf';
hhsuite on Biowulf


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



The HH-suite is an open-source software package for sensitive protein
sequence searching based on the pairwise alignment of hidden Markov models
(HMMs).



### References:


* M. Remmert, A. Biegert, A. Hauser, J. Söding. *HHblits: lightning-fast iterative protein sequence searching by HMM-HMM alignment.* Nature Methods 2011, 9:173-175. doi:[10.1038/nmeth.1818](https://doi.org/10.1038/nmeth.1818)


Documentation
* [hhsuite Main Site](https://github.com/soedinglab/hh-suite)
* [User Guide](https://github.com/soedinglab/hh-suite/blob/master/hhsuite-userguide.pdf)


Important Notes
* Module Name: hhsuite (see [the modules page](/apps/modules.html) for more information)
* hhsuite is a multithreaded application. Make sure to match the number of cpus requested with the number of threads.
* Environment variables set 
	+ HHLIB* Example files in /fdb/hhsuite/test-data* Reference data in /fdb/hhsuite/


hhblits and hhsearch need to do many random file access and read
operations. The central file systems (i.e. /data, /scratch, or /fdb) will
not perform well under this type of load. This means that running against
the database directly stored on /fdb will not be performant. In addition
it may tax the file system enough to also slow down other user's
computations. We therefore recommend to copy the database to be searched
to [lscratch](https://hpc.nih.gov/docs/userguide.html#local).
In particular, nodes with SSD storage should be used. This means that the
ideal usage pattern for large hhblits/hhsearch jobs is to allocate a node
(or nodes) exclusively, copy the database to lscratch, and then run all
computations on that node.



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive --constraint=ssd3200 --gres=lscratch:100 --cpus-per-task=10**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load hhsuite**
[+] Loading gcc  7.3.0  ... 
[+] Loading openmpi 3.0.0  for GCC 7.3.0 
[+] Loading hhsuite  3.3.0 
[user@cn3144 ~]$ **cd /lscratch/$SLURM\_JOB\_ID**
[user@cn3144 46116226]$ **cp -r /fdb/hhsuite/uniprot20\_2016\_02 .**
[user@cn3144 46116226]$ **hhblits \
 -cpu $SLURM\_CPUS\_PER\_TASK \
 -i /fdb/hhsuite/test-data/query.a3m \
 -d uniprot20\_2016\_02/uniprot20\_2016\_02 \
 -o query.hhr**
- 17:04:47.198 INFO: Searching 8290206 column state sequences.

- 17:04:47.723 INFO: /fdb/hhsuite/test-data/query.a3m is in A2M, A3M or FASTA format

- 17:04:47.794 INFO: Iteration 1

- 17:04:48.727 INFO: Prefiltering database

- 17:06:47.000 INFO: HMMs passed 1st prefilter (gapless profile-profile alignment)  : 302479
...
...
[user@cn3144 46116226]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. hhblits.sh). For example:



```

#!/bin/sh
set -e

module load hhsuite
cd /lscratch/$SLURM_JOB_ID
cp -r /fdb/hhsuite/uniprot20_2016_02 .

hhblits -i /fdb/hhsuite/test-data/query.seq \
  -d ./uniprot20_2016_02/uniprot20_2016_02 \
  -cpu $SLURM_CPUS_PER_TASK -o test.hhr \
  -oa3m test.a3m -n 6

hhmake -i test.a3m -o test.hhm
addss.pl test.hhm test_addss.hhm -hmm
cp test* /path/to/output/dir

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=6 --gres=lscratch:100 [--mem=#] hhblits.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. hhblits.swarm). For example:



```

cp -r /fdb/hhsuite/uniprot20_2016_02 /lscratch/$SLURM_JOB_ID \
&& hhblits -i sample1.seq -d /lscratch/$SLURM_JOB_ID/uniprot20_2016_02/uniprot20_2016_02 -cpu $SLURM_CPUS_PER_TASK -o sample1.hhr
cp -r /fdb/hhsuite/uniprot20_2016_02 /lscratch/$SLURM_JOB_ID \
&& hhblits -i sample2.seq -d /lscratch/$SLURM_JOB_ID/uniprot20_2016_02/uniprot20_2016_02 -cpu $SLURM_CPUS_PER_TASK -o sample2.hhr
cp -r /fdb/hhsuite/uniprot20_2016_02 /lscratch/$SLURM_JOB_ID \
&& hhblits -i sample3.seq -d /lscratch/$SLURM_JOB_ID/uniprot20_2016_02/uniprot20_2016_02 -cpu $SLURM_CPUS_PER_TASK -o sample3.hhr
cp -r /fdb/hhsuite/uniprot20_2016_02 /lscratch/$SLURM_JOB_ID \
&& hhblits -i sample4.seq -d /lscratch/$SLURM_JOB_ID/uniprot20_2016_02/uniprot20_2016_02 -cpu $SLURM_CPUS_PER_TASK -o sample4.hhr

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f hhblits.swarm [-g #] -t 6 --gres=lscratch:100 --module hhsuite
```

where


|  |  |  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --gres lscratch:# Number of gigabytes of local scratch space required for each process
 | --module hhsuite Loads the hhsuite module for each subjob in the swarm 
 | |
 | |
 | |
 | |








