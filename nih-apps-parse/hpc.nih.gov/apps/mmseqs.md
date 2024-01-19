

document.querySelector('title').textContent = 'mmseqs on Biowulf';
mmseqs on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
 |



From the Documentation:




> 
>  MMseqs2 (Many-against-Many sequence searching) is a software suite to search and cluster huge protein and nucleotide sequence sets. MMseqs2 is open source GPL-licensed software implemented in C++ for Linux, MacOS, and (as beta version, via cygwin) Windows. The software is designed to run on multiple cores and servers and exhibits very good scalability. MMseqs2 can run 10000 times faster than BLAST. At 100 times its speed it achieves almost the same sensitivity. It can perform profile searches with the same sensitivity as PSI-BLAST at over 400 times its speed.


### References:


* Steinegger M and Soeding J. *MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets.* Nature Biotechnology 2017, doi: 10.1038/nbt.3988. [PubMed](https://pubmed.ncbi.nlm.nih.gov/29035372/) | [Journal](https://www.nature.com/articles/nbt.3988)
* Steinegger M and Soeding J. *Clustering huge protein sequence sets in linear time.* Nature Communications 2018, doi: 10.1038/s41467-018-04964-5. [PubMed](https://pubmed.ncbi.nlm.nih.gov/29959318/) | [PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6026198/) | [Journal](https://www.nature.com/articles/s41467-018-04964-5)


Documentation
* [mmseqs on GitHub](https://github.com/soedinglab/MMseqs2)


Important Notes
* Module Name: mmseqs (see [the modules page](/apps/modules.html) for more information)
* mmseqs is multithreaded. Please match the number of threads with the number of allocated CPUs
* Example files in `$MMSEQS_TEST_DATA`



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) to follow the mmseqs metagenomic pathogen
detection tutorial. Note that for larger databases you would want more memory.



```

[user@biowulf]$ **sinteractive --gres=lscratch:10 -c2 --mem=6g**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144]$ **cd /lscratch/$SLURM\_JOB\_ID**
[user@cn3144]$ **module load mmseqs**
[user@cn3144]$ **cp ${MMSEQS\_TEST\_DATA:-none}/\* .**

```

 Create a taxonomically annotated sequence database from a subset of uniprot



```

[user@cn3144]$ **gunzip uniprot\_sprot\_2018\_03\_mapping.tsv.gz**
[user@cn3144]$ **mmseqs createdb uniprot\_sprot\_2018\_03.fasta.gz uniprot\_sprot**
createdb uniprot_sprot_2018_03.fasta.gz uniprot_sprot

MMseqs Version:         aabc78c298f35cbc7a4136206c1a83adaa68695f
Database type           0
Shuffle input database  true
Createdb mode           0
Write lookup file       1
Offset of numeric ids   0
Compressed              0
Verbosity               3

Converting sequences
[556915] 2s 838ms
Time for merging to uniprot_sprot_h: 0h 0m 0s 174ms
Time for merging to uniprot_sprot: 0h 0m 0s 358ms
Database type: Aminoacid
Time for processing: 0h 0m 4s 278ms

[user@cn3144]$ **mmseqs createtaxdb uniprot\_sprot tmp --tax-mapping-file uniprot\_sprot\_2018\_03\_mapping.tsv**
Create directory tmp
createtaxdb uniprot_sprot tmp --tax-mapping-file uniprot_sprot_2018_03_mapping.tsv

MMseqs Version:         aabc78c298f35cbc7a4136206c1a83adaa68695f
NCBI tax dump directory
Taxonomy mapping file   uniprot_sprot_2018_03_mapping.tsv
Taxonomy mapping mode   0
Taxonomy db mode        1
Threads                 72
Verbosity               3

Download taxdump.tar.gz
  % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                                 Dload  Upload   Total   Spent    Left  Speed
100 56.1M  100 56.1M    0     0   128M      0 --:--:-- --:--:-- --:--:--  129M
Loading nodes file ... Done, got 2442778 nodes
Loading merged file ... Done, added 68544 merged nodes.
Loading names file ... Done
Init RMQ ...Done


```

Similarity search to transfer taxonomic information to the unknown reads



```

[user@cn3144]$ **mmseqs createdb mystery\_reads.fasta reads**
createdb mystery_reads.fasta reads

MMseqs Version:         aabc78c298f35cbc7a4136206c1a83adaa68695f
Database type           0
Shuffle input database  true
Createdb mode           0
Write lookup file       1
Offset of numeric ids   0
Compressed              0
Verbosity               3

Converting sequences
[410667] 0s 417ms
Time for merging to reads_h: 0h 0m 0s 57ms
Time for merging to reads: 0h 0m 0s 114ms
Database type: Nucleotide
Time for processing: 0h 0m 1s 111ms

[user@cn3144]$ **mmseqs taxonomy reads uniprot\_sprot lca\_result tmp -s 2**
taxonomy reads uniprot_sprot lca_result tmp -s 2

MMseqs Version:                         aabc78c298f35cbc7a4136206c1a83adaa68695f
ORF filter                              1
ORF filter e-value                      100
ORF filter sensitivity                  2
LCA mode                                3
[...snip...]
[user@cn3144]$ **mmseqs createtsv reads lca\_result lca.tsv** # convert results to readable .tsv file
[user@cn3144]$ **mmseqs taxonomyreport uniprot\_sprot lca\_result report.txt** # create a summary report
[user@cn3144]$ **head report.txt**
1.7143  278     278     no rank 0       unclassified
98.2857 15939   8       no rank 1       root
88.6909 14383   254     no rank 131567    cellular organisms
79.9346 12963   335     superkingdom    2           Bacteria
73.8793 11981   607     phylum  1224          Proteobacteria
49.3741 8007    208     class   28216           Betaproteobacteria
46.7287 7578    452     order   80840             Burkholderiales
40.6857 6598    593     family  119060              Burkholderiaceae
31.4916 5107    893     genus   48736                 Ralstonia
16.0449 2602    1       species 305                     Ralstonia solanacearum

[user@cn3144]$ **mmseqs taxonomyreport uniprot\_sprot lca\_result report.html --report-mode 1** # krona plot

[user@cn3144]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. mmseqs.sh), which uses the input file 'mmseqs.in'. For example:



```

#!/bin/bash
module load mmseqs/2-13-45111-219-gaabc78c
cd /lscratch/$SLURM_JOB_ID
module load mmseqs
cp ${MMSEQS_TEST_DATA:-none}/* .
gunzip uniprot_sprot_2018_03_mapping.tsv.gz
mmseqs createdb uniprot_sprot_2018_03.fasta.gz uniprot_sprot
mmseqs createtaxdb uniprot_sprot tmp --tax-mapping-file uniprot_sprot_2018_03_mapping.tsv
mmseqs createdb mystery_reads.fasta reads
mmseqs taxonomy reads uniprot_sprot lca_result tmp -s 2
mmseqs createtsv reads lca_result lca.tsv
mmseqs taxonomyreport uniprot_sprot lca_result report.txt
mmseqs taxonomyreport uniprot_sprot lca_result report.html --report-mode 1

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=2 --mem=6g mmseqs.sh
```







