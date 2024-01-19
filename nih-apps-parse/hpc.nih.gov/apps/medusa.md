

document.querySelector('title').textContent = 'medusa on Biowulf';
medusa on Biowulf


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



A draft genome scaffolder that uses multiple reference genomes in a graph-based approach.




Features
* MEDUSA (Multi-Draft based Scaffolder), an algorithm for genome scaffolding. MEDUSA exploits information obtained from a set of (draft or closed) genomes from related organisms to determine the correct order and orientation of the contigs. MEDUSA formalizes the scaffolding problem by means of a combinatorial optimization formulation on graphs and implements an efficient constant factor approximation algorithm to solve it. In contrast to currently used scaffolders, it does not require either prior knowledge on the microrganisms dataset under analysis (e.g. their phylogenetic relationships) or the availability of paired end read libraries.




### References:


* Bosi E, Donati B, Galardini M, Brunetti S, Sagot MF, Lió P, Crescenzi P, Fani R, Fondi M. *MeDuSa: a multi-draft based scaffolder.* Bioinformatics. 2015 Aug 1;31(15):2443-51. doi: 10.1093/bioinformatics/btv171. Epub 2015 Mar 25. PMID: 25810435.
 [PubMed](https://pubmed.ncbi.nlm.nih.gov/25810435/) | 
 [Journal](https://academic.oup.com/bioinformatics/article/31/15/2443/188083)


Documentation
* medusa Main Site:[Main Site](https://github.com/combogenomics/medusa)


Important Notes
* Module Name: medusa (see [the modules page](/apps/modules.html) for more information)
 * Current medusa command lines could be run as
 
```

	medusa --help
```
* Environment variables set 
	+ $MEDUSAPATH* Example files in $MEDUSA\_TEST\_DATA



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load medusa**
[user@cn3144 ~]$ **cp -r $MEDUSA\_TEST\_DATA/\* .**
[user@cn3144 ~]$ **medusa --help**
medusa --help
Medusa version 1.6
usage: java -jar medusa.jar -i inputfile -v
available options:
 -d                                    OPTIONAL PARAMETER;The option *-d*
                                       allows for the estimation of the
                                       distance between pairs of contigs
                                       based on the reference genome(s):
                                       in this case the scaffolded contigs
                                       will be separated by a number of N
                                       characters equal to this estimate.
                                       The estimated distances are also
                                       saved in the
                                       \_distanceTable file.
 By default the scaffolded contigs
 are separated by 100 Ns
 -f <> OPTIONAL PARAMETER; The option \*-f\*
 is optional and indicates the path
 to the comparison drafts folder
 -gexf OPTIONAL PARAMETER;Conting network
 and path cover are given in gexf
 format.
 -h Print this help and exist.
 -i <> REQUIRED PARAMETER;The option \*-i\*
 indicates the name of the target
 genome file.
 -n50 <> OPTIONAL PARAMETER; The option
 \*-n50\* allows the calculation of
 the N50 statistic on a FASTA file.
 In this case the usage is the
 following: java -jar medusa.jar
 -n50 . All the
 other options will be ignored.
 -o <> OPTIONAL PARAMETER; The option \*-o\*
 indicates the name of output fasta
 file.
 -random <> OPTIONAL PARAMETER;The option
 \*-random\* is available (not
 required). This option allows the
 user to run a given number of
 cleaning rounds and keep the best
 solution. Since the variability is
 small 5 rounds are usually
 sufficient to find the best score.
 -scriptPath <> OPTIONAL PARAMETER; The folder
 containing the medusa scripts.
 Default value: medusa\_scripts
 -v RECOMMENDED PARAMETER; The option
 \*-v\* (recommended) print on console
 the information given by the
 package MUMmer. This option is
 strongly suggested to understand if
 MUMmer is not running properly.
 -w2 OPTIONAL PARAMETER;The option \*-w2\*
 is optional and allows for a
 sequence similarity based weighting
 scheme. Using a different weighting
 scheme may lead to better results.

[user@cn3144 ~]$ **medusa -f reference\_genomes/ -i Rhodobacter\_target.fna -v** 
INPUT FILE:Rhodobacter\_target.fna
------------------------------------------------------------------------------------------------------------------------
Running MUMmer...done.
------------------------------------------------------------------------------------------------------------------------
Building the network...done.
------------------------------------------------------------------------------------------------------------------------
Cleaning the network...done.
------------------------------------------------------------------------------------------------------------------------
Scaffolds File saved: Rhodobacter\_target.fnaScaffold.fasta
------------------------------------------------------------------------------------------------------------------------
Number of scaffolds: 78 (singletons = 32, multi-contig scaffold = 46)
from 564 initial fragments.
Total length of the jointed fragments: 4224838
Computing N50 on 78 sequences.
N50: 143991.0
----------------------
Summary File saved: Rhodobacter\_target.fna\_SUMMARY

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. medusa.sh). For example:



```

#!/bin/bash
set -e
module load medusa
medusa -f reference_genomes/ -i Rhodobacter_target.fna -v
```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=2 --mem=2g medusa.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. medusa.swarm). For example:



```

cd dir1;medusa -f reference_genomes/ -i 1_target.fna -v
cd dir2;medusa -f reference_genomes/ -i 2_target.fna -v 
cd dir3;medusa -f reference_genomes/ -i 3_target.fna -v
cd dir4;medusa -f reference_genomes/ -i 4_target.fna -v


```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f medusa.swarm [-g #] [-t #] --module medusa
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module medusa Loads the medusa module for each subjob in the swarm
 | |
 | |
 | |








