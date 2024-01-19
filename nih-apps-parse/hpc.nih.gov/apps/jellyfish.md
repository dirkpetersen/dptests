

document.querySelector('title').textContent = 'jellyfish on Biowulf';
jellyfish on Biowulf


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


 Jellyfish counts k-mers in fasta or fastq files (and sam/bam/cram starting
at version 2.2.7). k-mer counts are saved in a binary format that can be
queried or dumped to text based format. 


### References:


* G. Marcais and C. Kingsford. 
 *A fast, lock-free approach for efficient parallel counting of occurrences of k-mers*. 
 Bioinformatics 2011, 27:764-770.
 [PubMed](https://www.ncbi.nlm.nih.gov/pubmed/21217122) | 
 [PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3051319/) | 
 [Journal](https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/btr011)


Documentation
* [GitHub](https://github.com/gmarcais/Jellyfish)
* [Manual [pdf]](http://www.genome.umd.edu/docs/JellyfishUserGuide.pdf)
* [Home page](http://www.genome.umd.edu/jellyfish.html)


Important Notes
* Module Name: jellyfish (see [the modules page](/apps/modules.html) for more information)
* jellyfish is a multithreaded application. Please match the number of threads to the 
number of allocated CPUs
* Example files can be found in `$JELLYFISH_TEST_DATA`
* Some versions of jellyfish have a but when counting 7- or 8-mers.



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --cpus-per-task=6 --mem=10g**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load jellyfish**
[user@cn3144 ~]$ **cp -L $JELLYFISH\_TEST\_DATA/\* .**
[user@cn3144 ~]$ **ls -lh**
total 445M
-rw-r--r-- 1 user group 387M Feb 24 13:59 ENCFF001NGB.bam
-rw-r--r-- 1 user group  56M Feb 24 13:59 ERR458495.fastq.gz

[user@cn3144 ~]$ **# count 16-mers in fastq file (S. cerevisiae RNA-Seq data)**
[user@cn3144 ~]$ **# jellyfish does not natively read compressed data - uncompress on the fly**
[user@cn3144 ~]$ **jellyfish count -t $SLURM\_CPUS\_PER\_TASK \
 -m 16 -s 10M -C -o 16mer.jf <(zcat ERR458495.fastq.gz)**
[user@cn3144 ~]$ **ls -lh 16mer.jf**
-rw-r--r-- 1 user group 51M Feb 26 07:54 16mer.jf
[user@cn3144 ~]$ **jellyfish stats 16mer.jf**
Unique:    3319471
Distinct:  6625582
Total:     38369851
Max_count: 2284

[user@cn3144 ~]$ **jellyfish dump -L 2 -o 16mer.fa 16mer.jf**
[user@cn3144 ~]$ **head 16mer.fa**
>259
AAAAAAAAAAAAAAAA
>3
CAATTTAGCCTTTCGC
>2
CTCATCCATGTGAAAA
>2
AAGTCAGGCACAAATC
>2
CTTCATTTTGCCACCA

[user@cn3144 ~]$ **jellyfish query 16mer.jf AGCCAATTTGACTTCA**
AGCCAATTTGACTTCA 45
[user@cn3144 ~]$ **jellyfish histo -t $SLURM\_CPUS\_PER\_TASK 16mer.jf > hist**
[user@cn3144 ~]$ **head -5 hist**
1 3319471
2 1114139
3 594139
4 359875
5 238815

[user@cn3144 ~]$ **# count 16mers in alinged mouse data in bam format**
[user@cn3144 ~]$ **jellyfish count -t $SLURM\_CPUS\_PER\_TASK -m 16 -s 10M -C \
 -o mouse\_16mers.jf --sam ENCFF001NGB.bam**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. jellyfish.sh), which uses the input file 'jellyfish.in'. For example:



```

#!/bin/bash
module load jellyfish/2.2.7 || exit 1
jellyfish count -t $SLURM_CPUS_PER_TASK -m 16 -s 10M -C \
    -o mouse_16mers.jf --sam ENCFF001NGB.bam

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=6 --mem=10g jellyfish.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. jellyfish.swarm). For example:



```

jellyfish count -t $SLURM_CPUS_PER_TASK -m 12 -s 10M -o sample1_12.jf --sam sample1.bam
jellyfish count -t $SLURM_CPUS_PER_TASK -m 12 -s 10M -o sample2_12.jf --sam sample2.bam
jellyfish count -t $SLURM_CPUS_PER_TASK -m 12 -s 10M -o sample3_12.jf --sam sample3.bam

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f jellyfish.swarm -g 10 -t 4 --module jellyfish/2.2.7
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g #  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t #  Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module jellyfish  Loads the jellyfish module for each subjob in the swarm 
 | |
 | |
 | |








