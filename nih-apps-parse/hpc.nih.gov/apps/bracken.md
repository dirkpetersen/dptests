

document.querySelector('title').textContent = 'HemTools: a collection of NGS pipelines and bioinformatic analyses.';
**bracken: estimating species abundance in metagenomics data.**


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



Bracken is a companion program to Kraken 1, KrakenUniq, or Kraken 2 While Kraken classifies reads to multiple levels in the taxonomic tree, Bracken allows estimation of abundance at a single level using those classifications (e.g. Bracken can estimate abundance of species within a sample).



### References:


* Jennifer Lu, Florian P. Breitwieser, Peter Thielen and Steven L. Salzberg   

*Bracken: estimating species abundance in metagenomics data*   

 PeerJ Comput. Sci. 3:e104; DOI 10.7717/peerj-cs.104.
* Jennifer Lu, Natalia Rincon, Derrick E. Wood, Florian P. Breitwieser,
Christopher Pockrandt, Ben Langmead, Steven L. Salzberg and Martin Steinegger   

*Metagenome analysis using the Kraken software suite*    

Nature Protocols, v.17 (December 2022), 2815–2839


Documentation
* [Bracken Home page](https://ccb.jhu.edu/software/bracken/)


Important Notes
* Module Name: bracken (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set
	+ **BRACKEN\_HOME**  installation directory
	+ **BRACKEN\_BIN**  executable directory
	+ **BRACKEN\_DIR**  source code directory
	+ **BRACKEN\_DATA**  configuration files and models directory
	+ **BRACKEN\_BENCHMARKS**  benchmark datasets directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=64g --gres=gpu:p100:1,lscratch:300 -c16**
[user@cn3104 ~]$ **module load bracken**
[+] Loading kraken  2.1.2                          
[+] Loading bracken 2.8
[user@cn3104 ~]$ **bracken -h**
/usr/local/apps/bracken/2.8/bin/bracken: illegal option -- h
Usage: bracken -d MY_DB -i INPUT -o OUTPUT -w OUTREPORT -r READ_LEN -l LEVEL -t THRESHOLD
  MY_DB          location of Kraken database
  INPUT          Kraken REPORT file to use for abundance estimation
  OUTPUT         file name for Bracken default output
  OUTREPORT      New Kraken REPORT output file with Bracken read estimates
  READ_LEN       read length to get all classifications for (default: 100)
  LEVEL          level to estimate abundance at [options: D,P,C,O,F,G,S,S1,etc] (default: S)
  THRESHOLD      number of reads required PRIOR to abundance estimation to perform reestimation (default: 0)
[user@cn3104 ~]$ **bracken-build -h**
/usr/local/apps/bracken/2.8/bin/bracken-build: illegal option -- h
Usage: bracken_build -k KMER_LEN -l READ_LEN -d MY_DB -x K_INSTALLATION -t THREADS
  KMER_LEN       kmer length used to build the kraken database (default: 35)
  THREADS        the number of threads to use when running kraken classification and the bracken scripts
  READ_LEN       read length to get all classifications for (default: 100)
  MY_DB          location of Kraken database
  K_INSTALLATION location of the installed kraken/kraken-build scripts (default assumes scripts can be run from the user path)
[user@cn3104 ~]$ **cd /lscratch/$SLURM\_JOB\_ID**
[user@cn3104 ~]$ **mkdir fa && ca fa**
[user@cn3104 ~]$ **wget --quiet --input-file $KRAKEN\_TEST\_DATA/genomes\_for\_custom\_db.urls**
[user@cn3104 ~]$ **ls**
GCA_000006745.1_ASM674v1_genomic.fna.gz
GCA_000006885.1_ASM688v1_genomic.fna.gz
GCA_000007045.1_ASM704v1_genomic.fna.gz
GCA_000007825.1_ASM782v1_genomic.fna.gz
GCA_000008005.1_ASM800v1_genomic.fna.gz
GCA_000009005.1_ASM900v1_genomic.fna.gz
GCA_000009585.1_ASM958v1_genomic.fna.gz
...
[user@cn3104 ~]$ **gunzip \*.gz**
[user@cn3104 ~]$ **ls**
GCA_000006745.1_ASM674v1_genomic.fna
GCA_000006885.1_ASM688v1_genomic.fna
GCA_000007045.1_ASM704v1_genomic.fna
GCA_000007825.1_ASM782v1_genomic.fna
GCA_000008005.1_ASM800v1_genomic.fna
GCA_000009005.1_ASM900v1_genomic.fna
...
[user@cn3104 ~]$ **kraken2-build --download-taxonomy --db custom\_db -t 16**
Downloading nucleotide gb accession to taxon map... done.
Downloading nucleotide wgs accession to taxon map... done.
Downloaded accession to taxon map(s)
Downloading taxonomy tree data... done.
Uncompressing taxonomy data...

```

A folder custom\_db with subfolder taxonomy has been be created.

```

[user@cn3104 ~]$ **ls -t**
custom_db
GCA_002000745.1_ASM200074v1_genomic.fna
GCA_001641045.1_ASM164104v1_genomic.fna
GCA_001518875.1_ASM151887v1_genomic.fna
GCA_001518775.1_ASM151877v1_genomic.fna
GCA_001975045.1_ASM197504v1_genomic.fna
...
[user@cn3104 ~]$ **cd ..**
[user@cn3104 ~]$ **for f in fa/\*.fna; do
 kraken2-build --add-to-library $f --db custom\_db;
done**
Masking low-complexity regions of new file... done.
Added "fa/GCA_000006745.1_ASM674v1_genomic.fna" to library (custom_db)
Masking low-complexity regions of new file... done.
Added "fa/GCA_000006885.1_ASM688v1_genomic.fna" to library (custom_db)
Masking low-complexity regions of new file... done.
...

```

The folder custom\_db with subfolder library has been be created..   


```

[user@cn3104 ~]$ **cd custom\_db**
[user@cn3104 ~]$ **ln -s ../fa/custom\_db/taxonomy**
[user@cn3104 ~]$ **cd ..**
[user@cn3104 ~]$ **kraken2-build --download-library bacteria --db custom\_db -t 16**
Step 1/2: Performing rsync file transfer of requested files

Rsync file transfer complete.
Step 2/2: Assigning taxonomic IDs to sequences
...
Processed 38595 projects (91046 sequences, 161.35 Gbp)... done.
All files processed, cleaning up extra sequence files... done, library complete.
Masking low-complexity regions of downloaded library...

[user@cn3104 ~]$ **bracken-build -d custom\_db -t 10 -k 50 -l 500 -x ${KRAKEN\_DB} -t 16**

```

Exit the application:   


```

[user@cn3104 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```





