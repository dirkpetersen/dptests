

document.querySelector('title').textContent = 'HMMER on Biowulf';
HMMER on Biowulf


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



[![hammer_sm](/images/hammer_sm.gif)](http://hmmer.janelia.org)
### Profile hidden Markov models for biological sequence analysis


Profile hidden Markov models (profile HMMs) can be used to do sensitive
database searching using statistical descriptions of a sequence family's
consensus. HMMER uses profile HMMs, and can be useful in situations like:
* if you are working with an evolutionarily diverse protein family, a BLAST
search with any individual sequence may not find the rest of the sequences in
the family.
* the top hits in a BLAST search are hypothetical sequences from genome
projects.
* your protein consists of several domains which are of different types.





### References:


HMMER (pronounced 'hammer', as in a more precise mining tool than BLAST) was
developed by Sean Eddy at Washington University in St. Louis. 


Documentation
* The HMMER website
is [hmmer.janelia.org](http://hmmer.janelia.org/).
* [HMMER User Guide](http://hmmer.janelia.org/#documentation)
(PDF)


Important Notes
* Module Name: hmmer (see [the modules page](/apps/modules.html) for more information)
* HMMER is a cpu-intensive program and is parallelized using threads, so
that each instance of hmmsearch or the other search programs can use all the cpus allocated on a
node.
* An MPI version is also available, but given that most Biowulf nodes have 56 CPUs, only very large jobs may benefit from using the MPI version.
* environment variables set * Example files in /usr/local/apps/hmmer/tutorial* Reference data in /fdb/fastadb/



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --cpus-per-task=16 --mem=10g**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load hmmer**

[user@cn3144 ~]$ **hmmsearch --cpu $SLURM\_CPUS\_PER\_TASK /usr/local/apps/hmmer/tutorial/globins4.hmm /fdb/fastadb/nr.fas**
# hmmsearch :: search profile(s) against a sequence database
# HMMER 3.1b2 (February 2015); http://hmmer.org/
# Copyright (C) 2015 Howard Hughes Medical Institute.
# Freely distributed under the GNU General Public License (GPLv3).
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# query HMM file:                  /usr/local/apps/hmmer/tutorial/globins4.hmm
# target sequence database:        /fdb/fastadb/nr.fas
# number of worker threads:        16
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Query:       globins4  [M=149]
[...]
  Alignments for each domain:
  == domain 1  score: 18.2 bits;  conditional E-value: 0.0014
        globins4  41 qefFekFkdLstedelkksadvkkHgkkvldAlsdalakld..ekleaklkdLselHakklkvdpkyfkllsevlvdvlaarlpkeftadvqaal 133
                     q++F++  +L+   ++   +     g+ + +A+++  +++d  + l ++++ ++++H ++++++ ++++++++ l+++l +  +  ft dv+ a
  WP_087017392.1  30 QRMFDHNPELKDIFNMSH-QRTGRQGVALFEAVAAYAKNIDnlGALTTAVERIAHKH-TSFNIQAEHYQIVGHHLIETLRELASDAFTKDVEEAW 122
                     677777777873333333.3345679999***********87889999*********.58*******************************9886 PP

        globins4 134 e 134

  WP_087017392.1 123 T 123
                     5 PP



Internal pipeline statistics summary:
-------------------------------------
Query model(s):                              1  (149 nodes)
Target sequences:                    145198384  (53182534605 residues searched)
Passed MSV filter:                   4463369  (0.0307398); expected 2903967.7 (0.02)
Passed bias filter:                  3672241  (0.0252912); expected 2903967.7 (0.02)
Passed Vit filter:                    255552  (0.00176002); expected 145198.4 (0.001)
Passed Fwd filter:                     17247  (0.000118782); expected 1452.0 (1e-05)
Initial search space (Z):          145198384  [actual number of targets]
Domain search space  (domZ):           10862  [number of targets reported over threshold]
# CPU time: 1306.00u 76.30s 00:23:02.30 Elapsed: 00:11:03.54
# Mc/sec: 11942.31

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. HMMER.sh). The following example uses the query sequence set globins4.hmm in the example directory, run against the NCBI nr protein database.


```

#!/bin/bash
set -e
module load hmmer
hmmsearch --cpu $SLURM_CPUS_PER_TASK /usr/local/apps/hmmer/tutorial/globins4.hmm /fdb/fastadb/nr.fas

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=16 --mem=10g HMMER.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. HMMER.swarm). For example:



```


hmmsearch --cpu $SLURM_CPUS_PER_TASK  file1.fas  /fdb/fastadb/nr.fas
hmmsearch --cpu $SLURM_CPUS_PER_TASK  file2.fas  /fdb/fastadb/nr.fas
hmmsearch --cpu $SLURM_CPUS_PER_TASK  file3.fas  /fdb/fastadb/nr.fas
[...]

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f HMMER.swarm -g 10 -t 32 --module hmmer
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module hmmer Loads the HMMER module for each subjob in the swarm 
 | |
 | |
 | |










