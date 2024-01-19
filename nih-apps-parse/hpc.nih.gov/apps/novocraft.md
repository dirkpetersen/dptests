

document.querySelector('title').textContent = 'Novocraft on Biowulf';
Novocraft on Biowulf
NovoSplice beta is now available for testing.




|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
[Swarm of jobs](#swarm) 
[MPI jobs](#mpi) 
[Novoindex](#index) 
[Benchmarks](#bench) 
 |


 Novoalign is an aligner for single-ended and paired-end reads from the Illumina 
 Genome Analyser. Novoalign finds global optimum alignments using full Needleman-Wunsch 
 algorithm with affine gap penalties whilst performing at the same or better 
 speed than aligners that are limited to two mismatches and no insertions 
 or deletions.


Novoalign indexes for some common genome assemblies such as hg18 and hg19 
 are available in /fdb/novoalign. If there are other genomes you 
 want indexed, please email staff@hpc.nih.gov


Documentation
* [Novocraft homepage, see tab for documentation](http://www.novocraft.com/)



Important Notes
* Module Name: novocraft (see [the 
 modules page](/apps/modules.html) for more information)
* Multithreaded/MPI
* Reference data in /fdb/novoalign/





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

[user@cn3144 ~]$ **sinteractive --cpus-per-task=4 --mem=10g**
[user@cn3144 ~]$ **novoalign -c $SLURM\_CPUS\_PER\_TASK -d celegans -f sim1.fastq sim1r.fastq -o SAM > out1.sam**
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```




Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit). Create a batch input file (e.g. novo.sh). For example:

 
```

#!/bin/bash
set -e

module load novocraft

# cd to the appropriate directory
cd /data/$USER/mydir

# generate an index file named 'celegans' for the sequence file elegans.dna.fa
novoindex celegans elegans.dna.fa

# align the reads in file s_1_sequence.txt against the indexed genome of C.Elegans.
novoalign -c $SLURM_CPUS_PER_TASK -f s_1_sequence.txt -d celegans -o SAM > out.sam
```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.


$ sbatch --cpus-per-task=4 --mem=10g myscript 
 The number assigned to '--cpus-per-task' will be passed to the $SLURM\_CPUS\_PER\_TASK 
 in the script automatically. User can adjust memory requirement based on 
 needs using --mem as in the example. 



```
 
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. novo.swarm). For example:



```

cd /data/$USER/novo1; novoalign -c $SLURM_CPUS_PER_TASK -d celegans -f sim1.fastq sim1r.fastq -o SAM > out1.sam
cd /data/$USER/novo2; novoalign -c $SLURM_CPUS_PER_TASK -d celegans -f sim2.fastq sim2r.fastq -o SAM > out2.sam
cd /data/$USER/novo3; novoalign -c $SLURM_CPUS_PER_TASK -d celegans -f sim3.fastq sim3r.fastq -o SAM > out3.sam

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f novo.swarm -g 12 -t 4 --module novocraft
```

where
 

|  |  |
| --- | --- |
| -g *#*  | Number of Gigabytes of memory required for each process (1 line in the swarm command file)
  |
| -t *#* | Number of threads/CPUs required for each process (1 line in the swarm command file).
  |
| --module  | Loads the module for each subjob in the swarm  |


Running a NovoalignMPI or NovoalignCSMPI batch job on Biowulf
1. Create a batch script along the lines of the one below:



```

#!/bin/bash 
# the file name is novoMPIScript

# load the latest version of novoalignMPI
module load novocraft
cd /data/$USER/mydir
Nodelist=`make-novo-nodelist`
mpiexec -envall -host $Nodelist -np $SLURM_NTASKS novoalignMPI -c $SLURM_CPUS_PER_TASK -d /fdb/novoalign/chr_all_mm10.nix -f infile1.fq infile2.fq > outputfile
```

3. submit job on the biowulf headnode:


```

biowulf $ sbatch --partition=multinode --nodes=2 --ntasks=4 --cpus-per-task=28 --constraint=x2695 mpi.script 
```

* 'make-novo-nodelist' script will gather node name information and assign to variable $Nodelist
* $SLURM\_NTASKS **and** $SLURM\_CPUS\_PER\_TASK will be automatically assigned by slurm based on the submission command.
* In this example, a master task was first started which then forked out 2 novoalignMPI tasks to each of the 2 nodes (3 novoalign tasks and one master task). Each novoalign task will use 28 threads. 
* Note that one node will run 2 tasks (28x2=56 threads) and the other node will run one master task and one novoalign task (only 28-30 threads since the master task is only listenning).
 * The --constraint flag ensures that the job will allocate the same type of nodes.
* The -envall flag is for passing variables to mpiexec.



Novoindex memory usage


Novoindex can use a lot of memory, so it is worthwhile estimating the memory usage before submitting the job, to prevent nodes with overloaded memory. (Thanks to Colin Hercus of Novocraft for this information). 


 The memory used for a indexed genome is  

 N/2 + 4(k+1) + 4N/s   

 where N is the length of the reference genome, k the index k-mer length and s the indexing step size. Note that the second term must be converted to the same units as the first and third. 


 For example, for a 6GB reference sequence, with default values k=15 and s=2, the index size would be   

 6G/2 + 416 +4\*6G/2 = 3G + 4G + 12G = 20G 


 It might be better to set the options as -k=15 -s=3 and then have index of ~ 15G   

 or -k=14 -s=3 for an index size of 13G. 


 Changing k&s can have an effect on run time so it might be worth testing with a few values to find the best memory/run time trade off.



Novoalign & NovoalignMPI benchmark



In general, parallel jobs should scale to at least 70% efficiency for the sake of other Biowulf users. One user using twice the resources to squeeze out 10% more performance may be keeping other users from working at all. The efficiency of a parallel job can be calculated as follows, where *e* is efficiency, *n* is the number of processors running the simulation, *t1* is the performance time running on one node and *tn* is the performance time running on n nodes.



```
       *e = t1/(n\*tn)*  
```

For example if a job benchmarks at 20 minutes when running on one node and at 5 minutes when running on 8 nodes, we can figure the efficiency of scaling like this:



```
       e= 20/(8*5) = 50% 
```

50% is way below the 70% guideline; this job does not scale well out to 8 nodes and would therefore use too many resources to justify any benefit. Indeed it may be the case that adding nodes would slow the actual over-all performance (wall-time) of the job. This type of benchmarking should be done on all systems you intend to simulate before long-term runs.


To find the most appropriate number of nodes for a specific type of job, it is essential to run one's own benchmarks. 




 
 The following novoalignMPI benchmark was done using paired reads, each 27 gb, against hg19 index file, using either x2650 or x2695 nodes.


Since a master process is considered a task, it could take up to one node depend on how many cpus-per-task is set. For example, if 3 tasks and 56 cpus-per-task is requested, then 1 node which run the master process will pratically do nothing and only 2 nodes have novoalign jobs running. 


The following table has been colaborated so that only nodes running novoalign jobs are counted.




|  |  |  |  |
| --- | --- | --- | --- |
| No. of Nodes | x2650 | x2695 | Efficiency |
| 1 | 42m |  | 100% |
| 2 | 21m |  | 100% |
| 4 | 12m |  | 87.5% |
| 8 | 7m |  | 75% |
| 16 | 4.5m |  | 58% |
| 1 |  | 26m | 100% |
| 2 |  | 13m | 100% |
| 4 |  | 8.5m | 76% |
| 8 |  | 5.5m | 59% |
| 16 |  | 4.5m | 36% |


Based on the benchmark above, user should not use more than 4 nodes for x2695 and 8 nodes for x2650.


 


The following novoalign (non-MPI) benchmark was done using paired reads, each 27 gb, against hg19 index file, using either x2650 or x2695.




| CPUS | x2650 | x2695 | Efficiency |
| --- | --- | --- | --- |
| 1 | 643m |  | 100 |
| 2 | 343m |  | 94 |
| 4 | 224m |  | 72 |
| 8 |  147m  |  | 55 |
| 16 | 79m |  | 51 |
| 32 | 46m |  | 44 |
| 1 |  | 600m | 100 |
| 2 |  | 491m | 61 |
| 4 |  | 283m | 53 |
| 8 |  | 148m | 51 |
| 16 |  | 81m | 46 |
| 32 |  | 49m | 38 |


 






