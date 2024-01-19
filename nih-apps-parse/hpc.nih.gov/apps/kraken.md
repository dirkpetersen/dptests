

document.querySelector('title').textContent = 'kraken on Biowulf';
kraken on Biowulf


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


Kraken is a system for assigning taxonomic labels to short DNA sequences,
usually obtained through metagenomic studies. It uses exact matches of kmers
against a reference database.


### References:


* Derrick E Wood and Steven L Salzberg.*Kraken: ultrafast metagenomic
 sequence classification using exact alignments.* Genome Biology
 15(2014), R46.
 [Pubmed](http://www.ncbi.nlm.nih.gov/pubmed/24580807) | 
 [PMC](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4053813/) | 
 [Journal](http://www.genomebiology.com/2014/15/3/R46)



Do not run Kraken with the database in the default `/fdb` location or a shared `/data`
 directory. It *must* be run with the database in `/lscratch` or `/dev/shm`.



Documentation
* [Home page](http://ccb.jhu.edu/software/kraken/)
* [GitHub](https://github.com/DerrickWood/kraken)
* [Manual](http://ccb.jhu.edu/software/kraken/MANUAL.html)


Important Notes
* Module Name: kraken (see [the modules page](/apps/modules.html) for more information)
* kraken is a multithreaded application
* Example data can be found in `$KRAKEN_TEST_DATA`
* Kraken databases can be found in **`/fdb/kraken/`** which is
 also the value of `$KRAKEN_DB_PATH`. The standard database
 contains genomes of archaea, bacteria, and viruses from NCBI.
* kraken has undergone a number of changes with version 2 (default) including a different database format.
 kraken2-compatible databases in `/fdb/kraken` end in `_kraken2`.
* To run Kraken version 2 efficiently, 
 you should allocate at least as much memory as the size of the database and copy the database to
 `/lscratch`. Shrunken databases with
 subsets of kmers are provided to accomodate running with smaller memory requirements. 
 The current standard Kraken 2 database is about 60GB in size.
* Classification should **not** make direct use
 of the databases in `/fdb/kraken` since this can cause some severe file system strain.
 Database should either be copied to memory or lscratch.



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. In
the following sample session we will use kraken to build a custom database and run a classification
with the custom data base. Shown below is the same example for kraken and kraken2. kraken2 differs from
kraken 1 in several important ways.



 $(document).ready(function() {
 $('ul.tabs li').unbind('click');
 $('ul.tabs li').click(function(){
 var tab\_id = $(this).attr('data-tab');

 $(this).parent().children().removeClass('current');
 $(this).parent().children().each(function () {
 $("#"+$(this).attr('data-tab')).removeClass('current');
 });

 $(this).addClass('current');
 $("#"+tab\_id).addClass('current');
 });
 });


* Version 1
* Version 2




```

[user@biowulf]$ **sinteractive --mem=24g --cpus-per-task=4 --gres=lscratch:100**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load kraken/1.1**

###
### build a custom kraken database with some 400+ genomes
###
[user@cn3144 ~]$ **cd /lscratch/$SLURM\_JOB\_ID**
[user@cn3144 ~]$ **mkdir fa**
[user@cn3144 ~]$ **cd fa**
[user@cn3144 ~]$ **wget --quiet --input-file $KRAKEN\_TEST\_DATA/genomes\_for\_custom\_db.urls**
[user@cn3144 ~]$ **gunzip \*.gz**
[user@cn3144 ~]$ **ls -lh | head -5**
total 1.3G
-rw-r--r-- 1 user group 3.9M May 22  2016 GCA_000006745.1_ASM674v1_genomic.fna
-rw-r--r-- 1 user group 2.1M May 26  2016 GCA_000006885.1_ASM688v1_genomic.fna
-rw-r--r-- 1 user group 2.0M May 27  2016 GCA_000007045.1_ASM704v1_genomic.fna
-rw-r--r-- 1 user group 5.3M Aug  4  2014 GCA_000007825.1_ASM782v1_genomic.fna

[user@cn3144 ~]$ **cd ..**
[user@cn3144 ~]$ **kraken-build --download-taxonomy --db custom\_db**
[user@cn3144 ~]$ **for f in fa/\*.fna; do
 kraken-build --add-to-library $f --db custom\_db
 done**
Added "fa/GCA_000006745.1_ASM674v1_genomic.fna" to library (custom_db)
Added "fa/GCA_000006885.1_ASM688v1_genomic.fna" to library (custom_db)
[...snip...]
Added "fa/GCA_900129335.1_PP_HGAG_QV30_2SC_T4_genomic.fna" to library (custom_db)
Added "fa/GCA_900185995.1_BK1071_genomic.fna" to library (custom_db)

[user@cn3144 ~]$ **kraken-build --build --threads $SLURM\_CPUS\_PER\_TASK --db custom\_db**
Kraken build set to minimize disk writes.
Creating k-mer set (step 1 of 6)...
Found jellyfish v1.1.11
Hash size not specified, using '253231360'
K-mer set created. [6m11.633s]
Skipping step 2, no database reduction requested.
Sorting k-mer set (step 3 of 6)...
K-mer set sorted. [2m20.480s]
Skipping step 4, GI number to seqID map now obsolete.
Creating seqID to taxID map (step 5 of 6)...
721 sequences mapped to taxa. [3m48.166s]
Setting LCAs in database (step 6 of 6)...
Finished processing 721 sequences
Database LCAs set. [19m6.065s]
Database construction complete. [Total: 31m26.381s]

[user@cn3144 ~]$ **kraken-build --clean --db custom\_db**
[user@cn3144 ~]$ **du -sh custom\_db**
10G     custom_db/
[user@cn3144 ~]$ **tree custom\_db**
custom_db/
|-- [user    14K]  accmap_file.tmp
|-- [user   8.0G]  database.idx
|-- [user   1.7G]  database.kdb
`-- [user   4.0K]  taxonomy
    |-- [user   142M]  names.dmp
    `-- [user   109M]  nodes.dmp

[user@cn3144 ~]$ **cp -r custom\_db /data/$USER/kraken\_databases**

###
### Classification
###
[user@cn3144 ~]$ **cp -r $KRAKEN\_TEST\_DATA/accuracy .**
[user@cn3144 ~]$ **cd accuracy**
[user@cn3144 ~]$ **ls -lh**
total 7.5M
-rw-rw-r-- 1 user group 2.2K Aug 14 09:55 ACCURACY.README
-rw-rw-r-- 1 user group 1.2M Aug 14 09:55 HiSeq_accuracy.fa
-rw-rw-r-- 1 user group 1.8M Aug 14 09:55 MiSeq_accuracy.fa
-rw-rw-r-- 1 user group 4.6M Aug 14 09:55 simBA5_accuracy.fa
-rw-rw-r-- 1 user group  396 Aug 14 09:55 taxonomyDocSum.xslt

```

Now we will use the in memory temp file system `/dev/shm`
to store the database for classification. Note that unlike `lscratch`,
`/dev/shm` is not automatically cleaned at exit and care has to
be taken to clean up explicitly at the end of a job.



```

[user@cn3144 ~]$ **cp -r custom\_db /dev/shm**
[user@cn3144 ~]$ **kraken --db /dev/shm/custom\_db --fasta-input \
 --output HiSeq.kraken HiSeq\_accuracy.fa**
10000 sequences (0.92 Mbp) processed in 0.422s (1421.8 Kseq/m, 131.15 Mbp/m).
  8140 sequences classified (81.40%)
  1860 sequences unclassified (18.60%)

[user@cn3144 ~]$ **head -5 HiSeq.kraken**
C       A_hydrophila_HiSeq.922  380703  101     644:22 380703:1 644:30 0:18
C       A_hydrophila_HiSeq.1263 644     101     1236:3 644:1 2:1 644:8 1236:8 644:50
C       A_hydrophila_HiSeq.2905 1448139 101     644:8 1448139:51 0:1 644:11
C       A_hydrophila_HiSeq.4866 1448139 101     644:61 1448139:10
C       A_hydrophila_HiSeq.7009 644     101     0:13 644:14 0:31 644:13

### CLEANUP
[user@cn3144 ~]$ **rm -rf /dev/shm/custom\_db**

```

The same approach of copying to `/dev/shm` can be taken with the databases
from `/fdb/kraken` as long as they fit. Note that `/dev/shm` as
a whole is limited to 50% of memory and per job by the memory allocation. For example,
using a reduced size 'minikraken' database obtained from the kraken maintainers:



```

[user@cn3144 ~]$ **cp -r /fdb/kraken/20180220\_standard\_8GB /dev/shm**
[user@cn3144 ~]$ **kraken --db /dev/shm/20180220\_standard\_8GB --fasta-input \
 --output HiSeq.kraken.2 HiSeq\_accuracy.fa**
10000 sequences (0.92 Mbp) processed in 0.602s (997.0 Kseq/m, 91.96 Mbp/m).
  7523 sequences classified (75.23%)
  2477 sequences unclassified (24.77%)

```

Column three is the predicted taxid. Let's get the name and rank of the
unique taxids from NCBI eutils and compare them to the predicted taxa:

```

[user@cn3144 ~]$ **cut -f3 HiSeq.kraken.2 | sort -u | grep -wv 0 > taxid.uniq**
[user@cn3144 ~]$ **esum="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"**
[user@cn3144 ~]$ **for taxid in $(cat taxid.uniq); do
 curl -s -G "${esum}?db=taxonomy&id=${taxid}" | xsltproc taxonomyDocSum.xslt -
 done 2> /dev/null > taxa.uniq**
[user@cn3144 ~]$ **join -1 3 -2 1 -t '|' \
 <(sort -k3,3 HiSeq.kraken.2 | tr '\t' '|') \
 <(paste -d '|' taxid.uniq taxa.uniq | sort -t'|' -k1,1) \
 | awk -F'|' '{printf("%s prediction: %s(%s)\n", $3, $6, $7)}' \
 | sort > HiSeq.kraken.pred**
[user@cn3144 ~]$ **tail -n -100 HiSeq.kraken.pred | head -n 5**
X_axonopodis_HiSeq.875312 prediction: Xanthomonas phaseoli pv. dieffenbachiae LMG 695,()
X_axonopodis_HiSeq.87911 prediction: Xanthomonas citri group,species group()
X_axonopodis_HiSeq.881987 prediction: Xanthomonas,genus()
X_axonopodis_HiSeq.882508 prediction: Xanthomonas phaseoli pv. dieffenbachiae LMG 695,()
X_axonopodis_HiSeq.882724 prediction: Xanthomonas,genus()

### CLEANUP
[user@cn3144 ~]$ **rm -rf /dev/shm/20180220\_standard\_8GB**
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```








```

[user@biowulf]$ **sinteractive --mem=36g --cpus-per-task=4 --gres=lscratch:100**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load kraken/2.1.2**

###
### build a custom kraken database with some 400+ genomes
###
[user@cn3144 ~]$ **cd /lscratch/$SLURM\_JOB\_ID**
[user@cn3144 ~]$ **mkdir fa**
[user@cn3144 ~]$ **cd fa**
[user@cn3144 ~]$ **wget --quiet --input-file $KRAKEN\_TEST\_DATA/genomes\_for\_custom\_db.urls**
[user@cn3144 ~]$ **gunzip \*.gz**
[user@cn3144 ~]$ **ls -lh | head -5**
total 1.3G
-rw-r--r-- 1 user group 3.9M May 22  2016 GCA_000006745.1_ASM674v1_genomic.fna
-rw-r--r-- 1 user group 2.1M May 26  2016 GCA_000006885.1_ASM688v1_genomic.fna
-rw-r--r-- 1 user group 2.0M May 27  2016 GCA_000007045.1_ASM704v1_genomic.fna
-rw-r--r-- 1 user group 5.3M Aug  4  2014 GCA_000007825.1_ASM782v1_genomic.fna

[user@cn3144 ~]$ **cd ..**
[user@cn3144 ~]$ **kraken2-build --download-taxonomy --db custom\_db**
[user@cn3144 ~]$ **for f in fa/\*.fna; do
 kraken2-build --add-to-library $f --db custom\_db
 done**
Added "fa/GCA_000006745.1_ASM674v1_genomic.fna" to library (custom_db)
Added "fa/GCA_000006885.1_ASM688v1_genomic.fna" to library (custom_db)
[...snip...]
Added "fa/GCA_900129335.1_PP_HGAG_QV30_2SC_T4_genomic.fna" to library (custom_db)
Added "fa/GCA_900185995.1_BK1071_genomic.fna" to library (custom_db)

# larger databases would need more threads to build
[user@cn3144 ~]$ **kraken2-build --build --threads 2 --db custom\_db**
Creating sequence ID to taxonomy ID map (step 1)...
Found 721/721 targets, searched through 203908352 accession IDs, search complete.
Sequence ID to taxonomy ID map complete. [18.874s]
Estimating required capacity (step 2)...
Estimated hash table requirement: 258445896 bytes
Capacity estimation complete. [22.086s]
Building database files (step 3)...
Taxonomy parsed and converted.
CHT created with 8 bits reserved for taxid.
Completed processing of 721 sequences, 1343692302 bp
Writing data to disk...  complete.
Database files completed. [1m7.904s]
Database construction complete. [Total: 1m48.906s]

[user@cn3144 ~]$ **kraken2-build --clean --db custom\_db**
Database disk usage: 34G
After cleaning, database uses 247M
[user@cn3144 ~]$ **tree custom\_db**
custom_db/
|-- [user   246M]  hash.k2d
|-- [user     48]  opts.k2d
|-- [user    16K]  taxo.k2d
`-- [user    16K]  unmapped.txt


[user@cn3144 ~]$ **cp -r custom\_db /data/$USER/kraken\_databases**

###
### Classification
###
[user@cn3144 ~]$ **cp -r $KRAKEN\_TEST\_DATA/accuracy .**
[user@cn3144 ~]$ **cd accuracy**
[user@cn3144 ~]$ **ls -lh**
total 7.5M
-rw-rw-r-- 1 user group 2.2K Aug 14 09:55 ACCURACY.README
-rw-rw-r-- 1 user group 1.2M Aug 14 09:55 HiSeq_accuracy.fa
-rw-rw-r-- 1 user group 1.8M Aug 14 09:55 MiSeq_accuracy.fa
-rw-rw-r-- 1 user group 4.6M Aug 14 09:55 simBA5_accuracy.fa
-rw-rw-r-- 1 user group  396 Aug 14 09:55 taxonomyDocSum.xslt

[user@cn3144 ~]$ **kraken2 --db ../custom\_db --output HiSeq.kraken HiSeq\_accuracy.fa**
Loading database information... done.
10000 sequences (0.92 Mbp) processed in 0.100s (5993.0 Kseq/m, 552.81 Mbp/m).
  8381 sequences classified (83.81%)
  1619 sequences unclassified (16.19%)


[user@cn3144 ~]$ **head -5 HiSeq.kraken**
C       A_hydrophila_HiSeq.922  644     101     644:56 0:3 644:4 0:4
C       A_hydrophila_HiSeq.1263 644     101     2:2 644:5 1236:3 2:4 1236:4 644:24 1236:5 644:1 1236:5 644:14
C       A_hydrophila_HiSeq.2905 1448139 101     644:7 1448139:28 644:5 1448139:5 644:4 1448139:10 644:8
C       A_hydrophila_HiSeq.4866 1448139 101     644:60 9
C       A_hydrophila_HiSeq.7009 644     101     0:5 644:19 0:3 644:9 0:18 644:13

[user@cn3144 ~]$ **cp -r /fdb/kraken/20210223\_standard\_16GB\_kraken2 ..**
[user@cn3144 ~]$ **kraken2 --db ../20210223\_standard\_16GB\_kraken2 \
 --output HiSeq.kraken.2 HiSeq\_accuracy.fa**
Loading database information... done.
10000 sequences (0.92 Mbp) processed in 0.148s (4062.8 Kseq/m, 374.76 Mbp/m).
  9033 sequences classified (90.33%)
  967 sequences unclassified (9.67%)

```

Column three is the predicted taxid. Let's get the name and rank of the
unique taxids from NCBI eutils and compare them to the predicted taxa:

```

[user@cn3144 ~]$ **cut -f3 HiSeq.kraken.2 | sort -u | grep -wv 0 > taxid.uniq**
[user@cn3144 ~]$ **esum="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"**
[user@cn3144 ~]$ **for taxid in $(cat taxid.uniq); do
 curl -s -G "${esum}?db=taxonomy&id=${taxid}" | xsltproc taxonomyDocSum.xslt -
 done 2> /dev/null > taxa.uniq**
[user@cn3144 ~]$ **join -1 3 -2 1 -t '|' \
 <(sort -k3,3 HiSeq.kraken.2 | tr '\t' '|') \
 <(paste -d '|' taxid.uniq taxa.uniq | sort -t'|' -k1,1) \
 | awk -F'|' '{printf("%s prediction: %s(%s)\n", $3, $6, $7)}' \
 | sort > HiSeq.kraken.pred**
[user@cn3144 ~]$ **tail -n -100 HiSeq.kraken.pred | head -n 5**
X_axonopodis_HiSeq.875312 prediction: Xanthomonas phaseoli pv. dieffenbachiae LMG 695,()
X_axonopodis_HiSeq.87911 prediction: Xanthomonas citri group,species group()
X_axonopodis_HiSeq.881987 prediction: Xanthomonas,genus()
X_axonopodis_HiSeq.882508 prediction: Xanthomonas phaseoli pv. dieffenbachiae LMG 695,()
X_axonopodis_HiSeq.882724 prediction: Xanthomonas,genus()

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```






Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. kraken.sh) which copies a reduced size standard database
to `/lscratch` (or `/dev/shm` for version 1) and runs several fastq samples against a local copy of the database.



* Version 1
* Version 2




```

#! /bin/bash
# USAGE: kraken.sh dbname infile [infile ...]

module load kraken/1.1
db=$1
dbname=$(basename $1)
shift
[[ ${1:-none} != "none" ]] || exit 1

# make sure the database in /dev/shm is cleaned up
trap 'rm -rf /dev/shm/${dbname}' EXIT

[[ -d ${db} ]] && cp -r $db /dev/shm || exit 1

for infile in "$@"; do
    [[ -f ${infile} ]] \
    && kraken --db /dev/shm/${dbname} --threads ${SLURM_CPUS_PER_TASK} \
        --fastq-input --gzip-compressed --output ${infile}.kraken
done

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```

sbatch --cpus-per-task=16 --mem=38g kraken.sh \
    /fdb/kraken/20171019_standard_32GB sample0[0-9].fastq.gz

```




```

#! /bin/bash
# USAGE: kraken2.sh dbname infile [infile ...]

module load kraken/2.1.2
db=$1
dbname=$(basename $1)
shift
[[ ${1:-none} != "none" ]] || exit 1
    
[[ -d ${db} ]] && cp -r $db /lscratch/$SLURM_JOB_ID || exit 1

for infile in "$@"; do
    [[ -f ${infile} ]] \
    && kraken2 --db /lscratch/$SLURM_JOB_ID/${dbname} --threads ${SLURM_CPUS_PER_TASK} \
        --output ${infile}.kraken ${infile}
done

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```

  sbatch --cpus-per-task=16 --mem=72g kraken2.sh \
      /fdb/kraken/20220803_standard_kraken2 sample0[0-9].fastq.gz
  
```


Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. kraken2.swarm) that uses the script above to process several
fastq files in each job



```

kraken2.sh /fdb/kraken/20220803_standard_kraken2 sample0[0-9].fastq.gz
kraken2.sh /fdb/kraken/20220803_standard_kraken2 sample1[0-9].fastq.gz
kraken2.sh /fdb/kraken/20220803_standard_kraken2 sample2[0-9].fastq.gz

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f kraken2.swarm -g 72 -t 16 --module kraken/2.1.2
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g #  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t #  Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module kraken  Loads the kraken module for each subjob in the swarm
 | |
 | |
 | |








