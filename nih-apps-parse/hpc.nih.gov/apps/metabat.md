

document.querySelector('title').textContent = 'Metabat on Biowulf';
Metabat on Biowulf


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



MetaBAT: A robust statistical framework for reconstructing genomes from metagenomic data



### References:


* Kang et al. MetaBAT, an efficient tool for accurately reconstructing single genomes from complex microbial communities. PeerJ, 2015.


Documentation
* [metabat Main Site](https://bitbucket.org/berkeleylab/metabat/src/master/)


Important Notes
* Module Name: metabat (see [the modules page](/apps/modules.html) for more information)
 * Multithreaded app
* Environment variables set 
	+ METABAT\_TESTDATA=/usr/local/apps/metabat/TEST\_DATA* Example files in /usr/local/apps/metabat/TEST\_DATA



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive -c 4 --mem 10g**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load metabat**

[user@cn3144 ~]$ **cp $METABAT\_TESTDATA/\* .**

[user@cn3144 ~]$ **runMetaBat.sh -t $SLURM\_CPUS\_PER\_TASK assembly.fa \*.bam**
Executing: 'jgi_summarize_bam_contig_depths --outputDepth assembly.fa.depth.txt --pairedContigs assembly.fa.paired.txt --minContigLength 1000 --minContigDepth 1  library1.sorted.bam library2.sorted.bam' at Fri Jun 21 15:02:54 EDT 2019
Output depth matrix to assembly.fa.depth.txt
Output pairedContigs lower triangle to assembly.fa.paired.txt
minContigLength: 1000
minContigDepth: 1
jgi_summarize_bam_contig_depths v2.13-29-g2e72973 2019-06-12T22:50:43
Output matrix to assembly.fa.depth.txt
0: Opening bam: library1.sorted.bam
1: Opening bam: library2.sorted.bam
Allocating pairedContigs matrix: 0 MB over 2 threads
Processing bam files
[...]
Closing last bam file
Finished
Finished jgi_summarize_bam_contig_depths at Fri Jun 21 15:09:03 EDT 2019
Creating depth file for metabat at Fri Jun 21 15:09:03 EDT 2019
Executing: 'metabat2  --inFile assembly.fa --outFile assembly.fa.metabat-bins/bin --abdFile assembly.fa.depth.txt' at Fri Jun 21 15:09:03 EDT 2019
MetaBAT 2 (v2.13-29-g2e72973) using minContig 2500, minCV 1.0, minCVSum 1.0, maxP 95%, minS 60, and maxEdges 200. 
40 bins (93197986 bases in total) formed.
Finished metabat2 at Fri Jun 21 15:09:04 EDT 2019

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. metabat.sh). For example:



```

#!/bin/bash
module load metabat
runMetaBat.sh -t $SLURM_CPUS_PER_TASK assembly.fa *.bam

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=4 --mem=20g metabat.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. metabat.swarm). For example:



```

cd dir1; runMetaBat.sh -t $SLURM_CPUS_PER_TASK assembly.fa *.bam
cd dir2; runMetaBat.sh -t $SLURM_CPUS_PER_TASK assembly.fa *.bam
cd dir3; runMetaBat.sh -t $SLURM_CPUS_PER_TASK assembly.fa *.bam

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f metabat.swarm -g 20 -t 8 --module metabat
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module metabat Loads the metabat module for each subjob in the swarm 
 | |
 | |
 | |








