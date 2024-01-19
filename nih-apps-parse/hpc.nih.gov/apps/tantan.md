

document.querySelector('title').textContent = 'tantan on Biowulf';
tantan on Biowulf


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




Tantan masks low complexity and short-period tandem repeats in nucleic acid and
protein sequences. It either lower-cases the selected regions or replaces them
with a configurable letter.



### References:


* M. C. Frith. *A new repeat-masking method enables specific detection 
 of homologous sequences*. Nucleic Acids Research 2011, 39:e23.
 [PubMed](https://www.ncbi.nlm.nih.gov/pubmed/21109538) | 
 [PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3045581/) | 
 [Journal](https://academic.oup.com/nar/article-lookup/doi/10.1093/nar/gkq1212)



Documentation
* [GitLab](https://gitlab.com/mcfrith/tantan)


Important Notes
* Module Name: tantan (see [the modules page](/apps/modules.html) for more information)
* Example files in `$TANTAN_TEST_DATA`



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

[user@cn3144]$ **module load tantan**
[user@cn3144]$ # mask repetitive regions with lower case letters (DNA)
[user@cn3144]$ **tantan $TANTAN\_TEST\_DATA/CM001971.fa > test.fa**
[user@cn3144]$ **head -n6 test.fa**
>CM001971.1 Homo sapiens mitochondrion, complete sequence, whole genome shotgun sequence
CAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGT
ATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGAGCC
GGAGCACCCTATGTCGCAGTATCTGTCTTTGATTCCTGCCTCATCCTATT
ATTTATCGCACCTACGTTCAATATTACAGACGAACATACTTACTAAAGCG
TGTTAATtaattaaTGCTTGTAGGACATAATAATAACAATTGAATGTCTG
       |-----| masked region
[user@cn3144]$ # same, but replace with 'N' insead
[user@cn3144]$ **tantan -x N $TANTAN\_TEST\_DATA/CM001971.fa > test.fa**
[user@cn3144]$ **head -n6 test.fa**
>CM001971.1 Homo sapiens mitochondrion, complete sequence, whole genome shotgun sequence
CAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGT
ATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGAGCC
GGAGCACCCTATGTCGCAGTATCTGTCTTTGATTCCTGCCTCATCCTATT
ATTTATCGCACCTACGTTCAATATTACAGACGAACATACTTACTAAAGCG
TGTTAATNNNNNNNTGCTTGTAGGACATAATAATAACAATTGAATGTCTG
       |-----| masked region

[user@cn3144]$ # and for proteins
[user@cn3144]$ **tantan -p $TANTAN\_TEST\_DATA/NP\_002102.fa > test.fa**
[user@cn3144]$ **head -n3 test.fa**
>NP_002102.4 huntingtin [Homo sapiens]
MATLEKLMKAFESLKSFQqqqqqqqqqqqqqqqqqqqqqqpppppppppp
pqlpqpppqaqpllpqpqppppppppppgpavaEEPLHRPKKELSATKKD

[user@cn3144]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. tantan.sh) similar to the following example:



```

#! /bin/bash

module load tantan/13 || exit 1
tantan -p human_genome.fa > human_genome_masked.fa

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=2 --mem=4g tantan.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. tantan.swarm). For example:



```

tantan -p seq1.fa > seq1_masked.fa
tantan -p seq2.fa > seq2_masked.fa
tantan -p seq3.fa > seq3_masked.fa

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f tantan.swarm -g 2 -t 1 -p 2 --module tantan/13
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g #  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t #  Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module tantan  Loads the tantan module for each subjob in the swarm 
 | |
 | |
 | |








