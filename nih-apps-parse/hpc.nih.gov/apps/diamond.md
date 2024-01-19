

document.querySelector('title').textContent = 'Diamond on Biowulf';
Diamond on Biowulf


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



DIAMOND is a sequence aligner for protein and translated DNA searches, designed for high performance analysis of big sequence data. The key features are:

 * Pairwise alignment of proteins and translated DNA at 500x-20,000x speed of BLAST.
* Frameshift alignments for long read analysis.
* Low resource requirements and suitable for running on standard desktops or laptops.
* Various output formats, including BLAST pairwise, tabular and XML, as well as taxonomic classification.





### References:


* [Buchfink, Benjamin, Chao Xie, and Daniel H. Huson. "Fast and sensitive protein alignment using DIAMOND." *Nature methods* 12.1 (2015): 59-60.](https://www.nature.com/articles/nmeth.3176)


Documentation
* [Diamond Main Site](https://github.com/bbuchfink/diamond)


Important Notes
* Module Name: diamond (see [the modules page](/apps/modules.html) for more information)
* Multithreaded app (use -p option)
* Overloads and uses all available CPUs by default
* Example files in /usr/local/apps/diamond/TEST\_DATA



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf ~]$ **sinteractive -c8 --mem=10g --gres=lscratch:10**
salloc.exe: Pending job allocation 12273309
salloc.exe: job 12273309 queued and waiting for resources
salloc.exe: job 12273309 has been allocated resources
salloc.exe: Granted job allocation 12273309
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn0885 are ready for job
srun: error: x11: no local DISPLAY defined, skipping
error: unable to open file /tmp/slurm-spank-x11.12273309.0
slurmstepd: error: x11: unable to read DISPLAY value

[user@cn0885 ~]$ **module load diamond**
[+] Loading diamond  2.0.8  on cn0885

[user@cn0885 ~]$ **cd /lscratch/$SLURM\_JOB\_ID**

[user@cn0885 12273309]$ **cp /usr/local/apps/diamond/TEST\_DATA/\* .**

[user@cn0885 12273309]$ **diamond makedb --in uniprot\_sprot.fasta.gz -d uniprot\_sprot -p $SLURM\_CPUS\_PER\_TASK**
diamond v2.0.8.146 (C) Max Planck Society for the Advancement of Science
Documentation, support and updates available at http://www.diamondsearch.org

#CPU threads: 8
Scoring parameters: (Matrix=BLOSUM62 Lambda=0.267 K=0.041 Penalties=11/1)
Database input file: uniprot_sprot.fasta.gz
Opening the database file...  [0s]
Loading sequences...  [3.297s]
Masking sequences...  [2.902s]
Writing sequences...  [0.348s]
Hashing sequences...  [0.118s]
Loading sequences...  [0s]
Writing trailer...  [0.066s]
Closing the input file...  [0.017s]
Closing the database file...  [0.013s]
Database hash = 7190f6d1af560ffacdb1351b89d36883
Processed 556568 sequences, 199530821 letters.
Total time = 6.765s

[user@cn0885 12273309]$ **diamond blastx -d uniprot\_sprot.dmnd -q reads.fna -p ${SLURM\_CPUS\_PER\_TASK} -o matches.m8**
diamond v2.0.8.146 (C) Max Planck Society for the Advancement of Science
Documentation, support and updates available at http://www.diamondsearch.org

#CPU threads: 8
Scoring parameters: (Matrix=BLOSUM62 Lambda=0.267 K=0.041 Penalties=11/1)
Temporary directory:
#Target sequences to report alignments for: 25
Opening the database...  [0.115s]
Database: uniprot_sprot.dmnd (type: Diamond database, sequences: 556568, letters: 199530821)
Block size = 2000000000
Opening the input file...  [0s]
Opening the output file...  [0s]
Loading query sequences...  [0s]
Masking queries...  [0.003s]
Building query seed set...  [0.001s]
The host system is detected to have 405 GB of RAM. It is recommended to increase the block size for better performance using these parameters : -b8 -c1
Algorithm: Query-indexed
Building query histograms...  [0s]
Allocating buffers...  [0s]
Loading reference sequences...  [0.565s]
Masking reference...  [2.883s]
Initializing temporary storage...  [0s]
Building reference histograms...  [0.154s]
Allocating buffers...  [0s]
Processing query block 1, reference block 1/1, shape 1/1.
Building reference seed array...  [0.154s]
Building query seed array...  [0s]
Computing hash join...  [0s]
Building seed filter...  [0s]
Searching alignments...  [0.003s]
Deallocating buffers...  [0s]
Clearing query masking...  [0s]
Computing alignments...  [0.005s]
Deallocating reference...  [0.036s]
Loading reference sequences...  [0s]
Deallocating buffers...  [0s]
Deallocating queries...  [0s]
Loading query sequences...  [0s]
Closing the input file...  [0s]
Closing the output file...  [0s]
Closing the database file...  [0.008s]
Deallocating taxonomy...  [0s]
Total time = 4.089s
Reported 25 pairwise alignments, 25 HSPs.
1 queries aligned.
The host system is detected to have 405 GB of RAM. It is recommended to increase the block size for better performance using these parameters : -b8 -c1

[user@cn0885 12273309]$ **exit**
exit
salloc.exe: Relinquishing job allocation 12273309
salloc.exe: Job allocation 12273309 has been revoked.

[user@biowulf ~]$ 

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. diamond.sh). For example:



```

#!/bin/bash
set -e
module load diamond
diamond blastx -d uniprot_sprot.dmnd -q reads.fna -p ${SLURM_CPUS_PER_TASK} -o matches.m8

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=8 --mem=10g diamond.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. diamond.swarm). For example:



```

diamond blastx -d db_name -q read1.fna -p ${SLURM_CPUS_PER_TASK} -o out1
diamond blastx -d db_name -q read2.fna -p ${SLURM_CPUS_PER_TASK} -o out2
diamond blastx -d db_name -q read3.fna -p ${SLURM_CPUS_PER_TASK} -o out3

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f diamond.swarm -g 10 -t 8 --module diamond
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module diamond Loads the diamond module for each subjob in the swarm 
 | |
 | |
 | |








