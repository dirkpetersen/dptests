

document.querySelector('title').textContent = 'metaphlan on Biowulf';
metaphlan on Biowulf


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




 From the metaphlan documentation:
 



> MetaPhlAn is a computational tool for profiling the composition of microbial 
>  communities (Bacteria, Archaea and Eukaryotes) from metagenomic shotgun sequencing data
>  (i.e. not 16S) with species-level. With the newly added StrainPhlAn module, it is now
>  possible to perform accurate strain-level microbial profiling.
> 
>  MetaPhlAn relies on ~1.1M unique clade-specific marker genes identified from ~100,000
>  reference genomes (~99,500 bacterial and archaeal and ~500 eukaryotic), allowing:
> 
>  * unambiguous taxonomic assignments;
> * accurate estimation of organismal relative abundance;
> * species-level resolution for bacteria, archaea, eukaryotes and viruses;
> * strain identification and tracking;
> * orders of magnitude speedups compared to existing methods;
> * metagenomic strain-level population genomics
> 
> 
> 


Metaphlan versions >3 include *strainphlan*, *phylophlan*, and *hclust2*.


### References:


* D. Tin Truong, E. A. Franzosa, T. L. Tickle, M. Scholz, G. Weingart, E. Pasolli, A. Tett, C. Huttenhower and N. Segata. *MetaPhlAn2 for enhanced metagenomic taxonomic profiling*. Nat Methods. 13:101 (2016).
 [PubMed](https://www.ncbi.nlm.nih.gov/pubmed/26418763) | 
 [Journal](https://www.nature.com/articles/nmeth.3589)
* D. Tin Truong, A. Tett, E. Pasolli, C. Huttenhower, and N. Segata. *Microbial strain-level population structure and genetic diversity from metagenomes*. Genome Res. 27:626-638 (2017)
 [PubMed](https://www.ncbi.nlm.nih.gov/pubmed/28167665) | 
 [PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5378180/) | 
 [Journal](https://genome.cshlp.org/content/27/4/626.long)



Documentation
* [metaphlan github](https://github.com/biobakery/MetaPhlAn)


Important Notes
* Module Name: metaphlan (see [the modules page](/apps/modules.html) for more information)
* metaphlan uses multithreading. Please match the number of CPUs to the number of
 allocated cpus.
* Example files in `$METAPHLAN_TEST_DATA`
* We no longer provide metaphlan reference data since metaphlan3 requires these to be writable. This is
 not feasible in an HPC environment. Instead, users must install their own database the first time they run
 metaphlan. For example,
 
```

    	metaphlan --install --bowtie2db /data/$USER/mp\_db
    
```

 Use the --bowtie2db option to point to where you want to install the database. Subsequent calls with the same
 option and argument will use the installed database.
* For the older version, metaphlan2, reference data is in /fdb/metaphlan/v20/



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive -c4 --mem=8g --gres=lscratch:10**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144]$ **cd /lscratch/$SLURM\_JOB\_ID**

[user@cn3144]$ **module load metaphlan**

[user@cn3144]$ **export METAPHLAN\_DB=/data/$USER/mp\_db**

[user@cn3144]$ **mkdir fasta**

[user@cn3144]$ **cp ${METAPHLAN\_TEST\_DATA:-none}/\*.fasta.gz fasta**

[user@cn3144]$ **ls -lh fasta**
-rw-r--r-- 1 user group 690K Nov  4 13:11 SRS014459-Stool.fasta.gz
-rw-r--r-- 1 user group 608K Nov  4 13:11 SRS014464-Anterior_nares.fasta.gz
-rw-r--r-- 1 user group 704K Nov  4 13:11 SRS014470-Tongue_dorsum.fasta.gz
-rw-r--r-- 1 user group 748K Nov  4 13:11 SRS014472-Buccal_mucosa.fasta.gz
-rw-r--r-- 1 user group 696K Nov  4 13:11 SRS014476-Supragingival_plaque.fasta.gz
-rw-r--r-- 1 user group 687K Nov  4 13:11 SRS014494-Posterior_fornix.fasta.gz

[user@cn3144]$ **for f in fasta/\*.gz
do name=$(basename $f .fasta.gz)
 metaphlan --nproc $SLURM\_CPUS\_PER\_TASK --bowtie2out ${name}.bt2.bz2 \
 --input\_type fasta --bowtie2db $METAPHLAN\_DB $f > ${name}\_profile.txt
done**

[user@cn3144]$ **merge\_metaphlan\_tables.py \*\_profile.txt > merged\_abundance\_table.txt**


[user@cn3144]$ **exit**
salloc.exe: Relinquishing job allocation 46116226

[user@biowulf]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. metaphlan.sh). For example:



```

#!/bin/bash
module load metaphlan || exit 1
METAPHLAN_DB=/data/$USER/mp_db
fastq=/path/to/sample.fastq.gz
name=$(dirname $fastq)/$(basename $fastq .fastq.gz)
metaphlan --nproc $SLURM_CPUS_PER_TASK --bowtie2out ${name}.bt2.bz2 \
          --input_type fastq --bowtie2db $METAPHLAN_DB $fastq > ${name}_profile.txt

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=6 --mem=12g metaphlan.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. metaphlan.swarm). For example:



```

metaphlan --nproc $SLURM_CPUS_PER_TASK --bowtie2out sample1.bt2.bz2 \
          --input_type fastq --bowtie2db $METAPHLAN_DB sample1.fastq.gz > sample1_profile.txt
metaphlan --nproc $SLURM_CPUS_PER_TASK --bowtie2out sample2.bt2.bz2 \
          --input_type fastq --bowtie2db $METAPHLAN_DB sample2.fastq.gz > sample2_profile.txt
metaphlan --nproc $SLURM_CPUS_PER_TASK --bowtie2out sample3.bt2.bz2 \
          --input_type fastq --bowtie2db $METAPHLAN_DB sample3.fastq.gz > sample3_profile.txt

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f metaphlan.swarm -g 12 -t 6 --module metaphlan
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g #  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t #  Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module metaphlan  Loads the metaphlan module for each subjob in the swarm 
 | |
 | |
 | |








