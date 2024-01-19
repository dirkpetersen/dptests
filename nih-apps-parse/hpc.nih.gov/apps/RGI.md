

document.querySelector('title').textContent = 'RGI: Resistance Gene Identifier ';
**RGI: a robust antimicrobial resistance gene predicting tool.**


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


  


RGI (Resistance Gene Identifier) is a robust antimicrobial resistance (AMR) gene predicting tool.
It is based on newly curated Comprehensive Antibiotic Research Database (CARD) and allows detection
detect AMR genes from thirteen genomes of Pseudomonas strains.



  

### References:


* Bharat Kwatra   

*In Silico Prediction and Comparison of Resistomes in Model Pseudomonas
Strains by Resistance Gene Identifier (RGI)*  

[bioRxiv](https://www.biorxiv.org/content/10.1101/2021.11.15.468576v1.full) preprint doi: https://doi.org/10.1101/2021.11.15.468576; posted November 15, 2021


Documentation
* [RGI on Github](https://github.com/arpcard/rgi)
* [CARD Database](https://card.mcmaster.ca)


Important Notes
* Module Name: RGI (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set
	+ **RGI\_HOME**  installation directory
	+ **RGI\_BIN**       bin directory
	+ **RGI\_SRC**       source directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive** 
[user@cn3107 ~]$ **module load RGI**
[+] Loading singularity  3.10.5  on cn4203
[+] Loading RGI  6.0.2

```

RGI is used to predict resistome(s) from protein or nucleotide data based on homology and SNP models.   
   

The tool uses data from the CARD database.   
   

**Usage:**  

- Select your input sequence (in FASTA format).   

- Select your input type (CONTIG or PROTEIN).   

- Select your alignment tool (DIAMOND or BLAST).   

- Specify if you want to include loose hits (YES or NO).   

- Specify if you want to remove temporary files (YES or NO).   

- Specify if you want to low quality predictions (YES or NO).   

- Select your data type (WGS, PLASMID, CHROMOSOME or NA).   

- Run the tool.   
   

**Output:**   

There are 2 different output files produced by rgi:   

- summary.txt: a tabular file of all detected resistance genes, one gene per line, and   

- report.json: a json version of summary.txt.

```

[user@cn3107 ~]$ **rgi -h** 
my_command= exec
usage: rgi  []
 commands are:
 ---------------------------------------------------------------------------------------
 Database
 ---------------------------------------------------------------------------------------
 auto\_load Automatically loads CARD database, annotations and k-mer database
 load Loads CARD database, annotations and k-mer database
 clean Removes BLAST databases and temporary files
 database Information on installed card database
 galaxy Galaxy project wrapper

 ---------------------------------------------------------------------------------------
 Genomic
 ---------------------------------------------------------------------------------------

 main Runs rgi application
 tab Creates a Tab-delimited from rgi results
 parser Creates categorical JSON files RGI wheel visualization
 heatmap Heatmap for multiple analysis

 ---------------------------------------------------------------------------------------
 Metagenomic
 ---------------------------------------------------------------------------------------
 bwt Align reads to CARD and in silico predicted allelic variants (beta)

 ---------------------------------------------------------------------------------------
 Baits validation
 ---------------------------------------------------------------------------------------
 tm Baits Melting Temperature

 ---------------------------------------------------------------------------------------
 Annotations
 ---------------------------------------------------------------------------------------
 card\_annotation Create fasta files with annotations from card.json
 wildcard\_annotation Create fasta files with annotations from variants
 baits\_annotation Create fasta files with annotations from baits (experimental)
 remove\_duplicates Removes duplicate sequences (experimental)

 ---------------------------------------------------------------------------------------
 Pathogen of origin
 ---------------------------------------------------------------------------------------

 kmer\_build Build AMR specific k-mers database used for pathogen of origin (beta)
 kmer\_query Query sequences against AMR k-mers database to predict pathogen of origin (beta)



Resistance Gene Identifier - 6.0.2

positional arguments:
 {main,tab,parser,load,auto\_load,clean,galaxy,database,bwt,tm,card\_annotation,wildcard\_annotation,baits\_annotation,remove\_duplicates,heatmap,kmer\_build,kmer\_query}
 Subcommand to run

optional arguments:
 -h, --help show this help message and exit

Use the Resistance Gene Identifier to predict resistome(s) from protein or
nucleotide data based on homology and SNP models. Check
https://card.mcmaster.ca/download for software and data updates. Receive email
notification of monthly CARD updates via the CARD Mailing List
(https://mailman.mcmaster.ca/mailman/listinfo/card-l)

[user@cn3107 ~]$ **wget https://github.com/arpcard/rgi/archive/refs/tags/6.0.2.tar.gz** 
[user@cn3107 ~]$ **tar -zxf 6.0.2.tar.gz && rm -f 6.0.2.tar.gz && cd rgi-6.0.2**
[user@cn3107 ~]$ **./test.sh** 
=================================== CHECK DEPENDENCIES VERSIONS ===================================
my\_command= exec
blastp: 2.9.0+
 Package: blast 2.9.0, build Jun 19 2020 17:01:47
my\_command= exec
/opt/conda/envs/rgi/bin/bowtie2-align-s version 2.4.2
64-bit
Built on default-bf91a638-95fa-4b77-97c5-abccd9855c3e
Mon Nov 2 17:44:37 UTC 2020
Compiler: gcc version 7.5.0 (crosstool-NG 1.24.0.131\_87df0e6\_dirty)
Options: -O3 -msse2 -funroll-loops -g3 -fvisibility-inlines-hidden -std=c++17 -fmessage-length=0 -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem /opt/conda/envs/rgi/include -fdebug-prefix-map=/opt/conda/conda-bld/bowtie2\_1604337828546/work=/usr/local/src/conda/bowtie2-2.4.2 -fdebug-prefix-map=/opt/conda/envs/rgi=/usr/local/src/conda-prefix -DPOPCNT\_CAPABILITY -DWITH\_TBB -std=c++11 -DNO\_SPINLOCK -DWITH\_QUEUELOCK=1
...
Sizeof {int, long, long long, void\*, size\_t, off\_t}: {4, 8, 8, 8, 8, 8}
my\_command= exec
diamond version 2.1.7
my\_command= exec
samtools 1.9
Using htslib 1.9
Copyright (C) 2018 Genome Research Ltd.
my\_command= exec

bamtools 2.5.1
Part of BamTools API and toolkit
Primary authors: Derek Barnett, Erik Garrison, Michael Stromberg
(c) 2009-2012 Marth Lab, Biology Dept., Boston College

my\_command= exec
bedtools v2.30.0
my\_command= exec
jellyfish 1.1.10
my\_command= exec
KMA-1.3.23
=================================== RGI EXECUTABLE LOCATION ===================================
/usr/local/apps/RGI/6.0.2/bin/rgi
...
=================================== DOWNLOAD CARD CANONICAL DATA ===================================--2023-06-06 20:38:45-- https://card.mcmaster.ca/latest/data
Resolving dtn10-e0 (dtn10-e0)... 10.1.200.192
Connecting to dtn10-e0 (dtn10-e0)|10.1.200.192|:3128... connected.
Proxy request sent, awaiting response... 200 OK
Length: 3896906 (3.7M) [application/x-bzip2]
Saving to: ‘card\_data.tar.bz2’

 0K .......... .......... .......... .......... .......... 1% 706K 5s
 50K .......... .......... .......... .......... .......... 2% 1.39M 4s
 100K .......... .......... .......... .......... .......... 3% 62.4M 3s
...
 3550K .......... .......... .......... .......... .......... 94% 242M 0s
 3600K .......... .......... .......... .......... .......... 95% 215M 0s
 3650K .......... .......... .......... .......... .......... 97% 240M 0s
 3700K .......... .......... .......... .......... .......... 98% 249M 0s
 3750K .......... .......... .......... .......... .......... 99% 249M 0s
 3800K ..... 100% 154M=0.3s

2023-06-06 20:38:46 (12.9 MB/s) - ‘card\_data.tar.bz2’ saved [3896906/3896906]

=================================== DOWNLOAD CARD VARIANTS DATA ===================================
--2023-06-06 20:38:47-- https://card.mcmaster.ca/download/6/prevalence-v4.0.0.tar.bz2
Resolving dtn10-e0 (dtn10-e0)... 10.1.200.192
Connecting to dtn10-e0 (dtn10-e0)|10.1.200.192|:3128... connected.
Proxy request sent, awaiting response... 200 OK
Length: 1537517791 (1.4G) [application/x-bzip2]
Saving to: ‘prevalence-v4.0.0.tar.bz2’

 0K .......... .......... .......... .......... .......... 0% 701K 35m41s
 50K .......... .......... .......... .......... .......... 0% 1.39M 26m40s
 100K .......... .......... .......... .......... .......... 0% 42.6M 17m58s
...
1501000K .......... .......... .......... .......... .......... 99% 227M 0s
1501050K .......... .......... .......... .......... .......... 99% 268M 0s
1501100K .......... .......... .......... .......... .......... 99% 260M 0s
1501150K .......... .......... .......... .......... .......... 99% 274M 0s
1501200K .......... .......... .......... .......... .......... 99% 237M 0s
1501250K .......... .......... .......... .......... .......... 99% 264M 0s
1501300K .......... .......... .......... .......... .......... 99% 258M 0s
1501350K .......... .......... .......... .......... .......... 99% 260M 0s
1501400K .......... .......... .......... .......... .......... 99% 233M 0s
1501450K .......... .......... .......... .. 100% 269M=16s

2023-06-06 20:39:04 (89.5 MB/s) - ‘prevalence-v4.0.0.tar.bz2’ saved [1537517791/1537517791]

=================================== CARD CANONICAL ANNOTATIONS ===================================
COMMAND rgi card\_annotation --input card\_data/card.json
my\_command= exec
=================================== VERSIONS ===================================
COMMAND data\_version: 3.2.7
COMMAND variants\_version: 4.0.0
=================================== CARD VARIANTS ANNOTATIONS ===================================
COMMAND rgi wildcard\_annotation --input\_directory card\_variants --version '4.0.0' --card\_json card\_data/card.json
my\_command= exec
=================================== CLEAN OLD DATABASES ===================================
COMMAND rgi clean --debug
my\_command= exec
INFO 2023-06-07 14:29:59,649 : Remove: /vf/users/user/RGI/rgi-6.0.2/app/\_db/16s\_rRNA.txt
INFO 2023-06-07 14:29:59,652 : Remove: /vf/users/user/RGI/rgi-6.0.2/app/\_db/23s\_rRNA.txt
...
INFO 2023-06-07 14:30:00,144 : Remove: /vf/users/user/RGI/rgi-6.0.2/app/\_data/variants
INFO 2023-06-07 14:30:00,365 : Cleaned directory: /vf/users/user/RGI/rgi-6.0.2/app/\_data/
INFO 2023-06-07 14:30:00,367 : Cleaned directory: /vf/users/user/RGI/rgi-6.0.2/app/\_db/
=================================== LOAD DATABASES ===================================
COMMAND rgi load --card\_json card\_data/card.json --card\_annotation card\_database\_v3.2.7.fasta --wildcard\_index card\_variants/index-for-model-sequences.txt --wildcard\_version '4.0.0' --wildcard\_annotation wildcard\_database\_v4.0.0.fasta --debug
...

```

End the interactive session:

```

[user@cn3107 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```





