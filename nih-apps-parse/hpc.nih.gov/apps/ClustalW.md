

document.querySelector('title').textContent = 'ClustalW on Biowulf';
ClustalW on Biowulf


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



Clustal W is a general purpose multiple alignment program for DNA or proteins. The sensitivity of the commonly used progressive multiple sequence alignment method has been greatly improved for the alignment of divergent protein sequences. It is designed to be run interactively, or to assign options via the command line.

Clustal Omega is a new development to the Clustal family, which offers a significant increase in scalability over previous versions, allowing hundreds of thousands of sequences to be aligned in only a few hours. It will also make use of multiple processors, where present. In addition, the quality of alignments is superior to previous versions, as measured by a range of popular benchmarks.

ClustalW is no longer being maintained or updated by its developers. ClustalO is actively being maintained. ClustalO can also run in multi-threaded mode, which may make your sequence alignments faster. 



Documentation
* [ClustalW website](http://www.clustal.org/clustal2/)
* [Clustal-Omega website](http://www.clustal.org/omega/)


Important Notes
* ClustalW Module Name: clustalw (see [the modules page](/apps/modules.html) for more information)
* Clustal-Omega Module Name: clustalo (see [the modules page](/apps/modules.html) for more information)
* ClustalW is single-threaded, while ClustalO can be multithreaded.
* Environment variables set



[A web interface to ClustalW](http://helixweb.nih.gov/multi-align/) and other multiple sequence alignment programs is available on 
our systems. It has been developed in-house, and allows the user to do a multiple sequence alignment using several different programs and
compare the results. 

Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive --cpus-per-task=8**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load clustalo**
[+] Loading clustalo 1.2.4  ...

[user@cn3144 ~]$ **clustalo -i globins630.fa -o clustalo.out --threads=$SLURM\_CPUS\_PER\_TASK**

[user@cn3144 ~]$ **module load clustalw**

[user@cn3144 ~]$ **clustalw**
**************************************************************
 ******** CLUSTAL W (1.83) Multiple Sequence Alignments  ********
 **************************************************************

     1. Sequence Input From Disc
     2. Multiple Alignments
     3. Profile / Structure Alignments
     4. Phylogenetic trees

     S. Execute a system command
     H. HELP
     X. EXIT (leave program)

Your choice: **1**

Sequences should all be in 1 file.

7 formats accepted: 
NBRF/PIR, EMBL/SwissProt, Pearson (Fasta), GDE, Clustal, GCG/MSF, RSF.

Enter the name of the sequence file: **seqs.inp**

Sequence format is Pearson
Sequences assumed to be PROTEIN

Sequence 1: chiins          110 aa
Sequence 2: xenins          110 aa
Sequence 3: humins          110 aa
Sequence 4: monins          110 aa
Sequence 5: dogins          110 aa
Sequence 6: hamins          110 aa
Sequence 7: bovins          110 aa
Sequence 8: guiins          110 aa

 **************************************************************
 ******** CLUSTAL W (1.83) Multiple Sequence Alignments  ********
 **************************************************************

     1. Sequence Input From Disc
     2. Multiple Alignments
     3. Profile / Structure Alignments
     4. Phylogenetic trees

     S. Execute a system command
     H. HELP
     X. EXIT (leave program)

Your choice: **2**
****** MULTIPLE ALIGNMENT MENU ******
    1.  Do complete multiple alignment now (Slow/Accurate)
    2.  Produce guide tree file only
    3.  Do alignment using old guide tree file

    4.  Toggle Slow/Fast pairwise alignments = SLOW

    5.  Pairwise alignment parameters
    6.  Multiple alignment parameters

    7.  Reset gaps before alignment? = OFF
    8.  Toggle screen display          = ON
    9.  Output format options

    S.  Execute a system command
    H.  HELP
    or press [RETURN] to go back to main menu

Your choice: **1**

Enter a name for the CLUSTAL output file  [seqs.aln]: **myseqs.aln**

Enter name for new GUIDE TREE           file   [seqs.dnd]: 

Start of Pairwise alignments
Aligning...
Sequences (1:2) Aligned. Score:  66
Sequences (1:3) Aligned. Score:  63
[...]
Sequences (7:8) Aligned. Score:  61
Guide tree        file created:   [seqs.dnd]
Start of Multiple Alignment
There are 7 groups
Aligning...
Group 1: Sequences:   2      Score:2138
Group 2: Sequences:   2      Score:2373
Group 3: Sequences:   4      Score:2157
Group 4: Sequences:   5      Score:2146
Group 5: Sequences:   6      Score:1971
Group 6: Sequences:   2      Score:1972
Group 7: Sequences:   8      Score:1739
Alignment Score 13222

Consensus length = 110
CLUSTAL-Alignment file created  [seqs.aln]

CLUSTAL W (1.83) multiple sequence alignment

dogins          MALWMRLLPLLALLALWAPAPTRAFVNQHLCGSHLVEALYLVCGERGFFYTPKARREVED
bovins          MALWTRLRPLLALLALWPPPPARAFVNQHLCGSHLVEALYLVCGERGFFYTPKARREVEG
humins          BALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAED
monins          BALWMRLLPLLALLALWGPDPVPAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAED
hamins          MTLWMRLLPLLTLLVLWEPNPAQAFVNQHLCGSHLVEALYLVCGERGFFYTPKSRRGVED
guiins          MALWMHLLTVLALLALWGPNTGQAFVSRHLCGSNLVETLYSVCQDDGFFYIPKDRRELED
chiins          BALWIRSLPLLALLVFSGPGTSYAAANQHLCGSHLVEALYLVCGERGFFYSPKARRDVEQ
xenins          BALWMQCLPLVLVLFFSTPNTE-ALVNQHLCGSHLVEALYLVCGDRGFFYYPKVKRDMEQ
                 :** :  .:: :* :  * .  * ..:*****:***:** ** : **** ** :*  * 

dogins          LQVRDVELAGAPGEGGLQPLALEGALQKRGIVEQCCTSICSLYQLENYCN
bovins          PQVGALELAGGPGAGG-----LEGPPQKRGIVEQCCASVCSLYQLENYCN
humins          LQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN
monins          PQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN
hamins          PQVAQLELGGGPGADDLQTLALEVAQQKRGIVDQCCTSICSLYQLENYCN
guiins          PQVEQTELGMGLGAGGLQPLALEMALQKRGIVDQCCTGTCTRHQLQSYCN
chiins          PLVSS---PLRGEAGVLPFQQEEYEKVKRGIVEQCCHNTCSLYQLENYCN
xenins          ALVSG---PQDNELDGMQLQPQEYQKMKRGIVEQCCHSTCSLFQLESYCN
                  *           .       *    *****:*** . *: .**:.***

Press [RETURN] to continue or  X  to stop: X

 **************************************************************
 ******** CLUSTAL W (1.83) Multiple Sequence Alignments  ********
 **************************************************************

     1. Sequence Input From Disc
     2. Multiple Alignments
     3. Profile / Structure Alignments
     4. Phylogenetic trees

     S. Execute a system command
     H. HELP
     X. EXIT (leave program)

Your choice: x

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. clustal.sh). For example, the following script runs ClustalW and ClustalO on the same input file.



```

#!/bin/bash

module load clustalw
clustalw -INFILE=test1.fa -ALIGN 

module load clustalo
clustalo -i test1.fa -o test1.out --threads=$SLURM_CPUS_PER_TASK

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=8 [--mem=#] clustalw.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. clustalw.swarm). For example:



```

clustalo -i test1.fa -o test1.out --threads=$SLURM_CPUS_PER_TASK
clustalo -i test2.fa -o test2.out --threads=$SLURM_CPUS_PER_TASK
clustalo -i test3.fa -o test3.out --threads=$SLURM_CPUS_PER_TASK
[...]

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f clustalw.swarm -t 8  [-g #] --module clustalw
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module clustalw Loads the clustalw module for each subjob in the swarm 
 | |
 | |
 | |










