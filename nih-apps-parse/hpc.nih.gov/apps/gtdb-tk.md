

document.querySelector('title').textContent = "gtdb-tk";
GTDB-TK on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
 |



From the GTDB-TK documentation




> 
>  GTDB-Tk is a software toolkit for assigning objective taxonomic
>  classifications to bacterial and archaeal genomes based on the Genome
>  Database Taxonomy ([GTDB](https://gtdb.ecogenomic.org/))). It is
>  designed to work with recent advances that allow hundreds or thousands of
>  metagenome-assembled genomes (MAGs) to be obtained directly from
>  environmental samples. It can also be applied to isolate and single-cell
>  genomes.
> 


### References:


* Chaumeil P. A. et al.
 *GTDB-Tk v2: memory friendly classification with the Genome Taxonomy Database*
 Bioinformatics 2022.
 [PubMed](https://pubmed.ncbi.nlm.nih.gov/36218463/) | 
 [PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/pmid/36218463/) | 
 [Journal](https://academic.oup.com/bioinformatics/article/38/23/5315/6758240)
* Chaumeil P. A. et al.
 *GTDB-Tk: A toolkit to classify genomes with the Genome Taxonomy Database*
 Bioinformatics 2019
 [PubMed](https://pubmed.ncbi.nlm.nih.gov/31730192/) | 
 [PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/pmid/31730192/) | 
 [Journal](https://academic.oup.com/bioinformatics/article/36/6/1925/5626182)


Documentation
* [Manual](https://ecogenomics.github.io/GTDBTk)
* [GitHub repository](https://github.com/Ecogenomics/GTDBTk)


Important Notes
* Module Name: gtdb-tk (see [the
 modules page](/apps/modules.html) for more information)
* GTDB-TK can use multiple CPUs. Please match allocated CPUs with `--cpu`. More
 CPUs may require more memory and the tool scales efficiently to no more than 16 CPUs.
* While 50GB were sufficient for the example below your data may require more memory. In
 one example with a large number of input genomes pplacer required ~420GB to complete.
 The documentation suggests that reducing the number of CPUs for pplacer may help but in
 our testing that did not appear to be the case.
* Please allocate lscratch and use `--tmpdir`
* The output folder should probably be in lscratch and then moved to /data at the end
* Example files in `$GTDBTK_TEST_DATA`
* Reference data in `/fdb/gtdb-tk`
* Environment variables: `$GTDBTK_DB`



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive --mem=50g --cpus-per-task=2 --gres=lscratch:50**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load gtdb-tk/2.3.0**

[user@cn3144 ~]$ **cp -r $GTDBTK\_TEST\_DATA genomes**
[user@cn3144 ~]$ **mkdir /lscratch/$SLURM\_JOB\_ID/tmp**
[user@cn3144 ~]$ **gtdbtk identify --genome\_dir genomes \
 --tmpdir /lscratch/$SLURM\_JOB\_ID/tmp \
 --out\_dir /lscratch/$SLURM\_JOB\_ID/identify \
 --extension gz \
 --cpus $SLURM\_CPUS\_PER\_TASK**
[2023-07-06 09:53:07] INFO: GTDB-Tk v2.3.0
[2023-07-06 09:53:07] INFO: gtdbtk identify --genome_dir genomes --out_dir identify --extension gz --cpus 2
[2023-07-06 09:53:07] INFO: Using GTDB-Tk reference data version r214: /refdata/
[2023-07-06 09:53:07] INFO: Identifying markers in 2 genomes with 2 threads.
[2023-07-06 09:53:07] TASK: Running Prodigal V2.6.3 to identify genes.
[2023-07-06 09:53:14] INFO: Completed 2 genomes in 7.21 seconds (3.61 seconds/genome).
[2023-07-06 09:53:14] TASK: Identifying TIGRFAM protein families.
[2023-07-06 09:53:19] INFO: Completed 2 genomes in 4.90 seconds (2.45 seconds/genome).
[2023-07-06 09:53:19] TASK: Identifying Pfam protein families.
[2023-07-06 09:53:19] INFO: Completed 2 genomes in 0.35 seconds (5.67 genomes/second).
[2023-07-06 09:53:19] INFO: Annotations done using HMMER 3.3.2 (Nov 2020).
[2023-07-06 09:53:19] TASK: Summarising identified marker genes.
[2023-07-06 09:53:19] INFO: Completed 2 genomes in 0.04 seconds (55.76 genomes/second).
[2023-07-06 09:53:19] INFO: Done.

[user@cn3144 ~]$ **gtdbtk align --identify\_dir /lscratch/$SLURM\_JOB\_ID/identify \
 --out\_dir /lscratch/$SLURM\_JOB\_ID/align \
 --tmpdir /lscratch/$SLURM\_JOB\_ID/tmp \
 --cpus $SLURM\_CPUS\_PER\_TASK**
[2023-07-06 09:54:01] INFO: GTDB-Tk v2.3.0
[2023-07-06 09:54:01] INFO: gtdbtk align --identify_dir identify --out_dir align --cpus 2
[2023-07-06 09:54:01] INFO: Using GTDB-Tk reference data version r214: /refdata/
[2023-07-06 09:54:01] INFO: Aligning markers in 2 genomes with 2 CPUs.
[2023-07-06 09:54:01] INFO: Processing 2 genomes identified as archaeal.
[2023-07-06 09:54:01] INFO: Read concatenated alignment for 4,416 GTDB genomes.
[2023-07-06 09:54:01] TASK: Generating concatenated alignment for each marker.
[2023-07-06 09:54:01] INFO: Completed 2 genomes in 0.02 seconds (115.54 genomes/second).
[2023-07-06 09:54:01] TASK: Aligning 52 identified markers using hmmalign 3.3.2 (Nov 2020).
[2023-07-06 09:54:02] INFO: Completed 52 markers in 0.32 seconds (162.74 markers/second).
[2023-07-06 09:54:02] TASK: Masking columns of archaeal multiple sequence alignment using canonical mask.
[2023-07-06 09:54:05] INFO: Completed 4,418 sequences in 3.77 seconds (1,171.18 sequences/second).
[2023-07-06 09:54:05] INFO: Masked archaeal alignment from 13,540 to 10,135 AAs.
[2023-07-06 09:54:05] INFO: 0 archaeal user genomes have amino acids in <10.0% of columns in filtered MSA.
[2023-07-06 09:54:05] INFO: Creating concatenated alignment for 4,418 archaeal GTDB and user genomes.
[2023-07-06 09:54:07] INFO: Creating concatenated alignment for 2 archaeal user genomes.
[2023-07-06 09:54:07] INFO: Done.

[user@cn3144 ~]$ **gtdbtk classify --genome\_dir genomes \
 --align\_dir /lscratch/$SLURM\_JOB\_ID/align \
 --out\_dir /lscratch/$SLURM\_JOB\_ID/classify \
 --tmpdir /lscratch/$SLURM\_JOB\_ID/tmp \
 -x gz \
 --cpus $SLURM\_CPUS\_PER\_TASK \
 --skip\_ani\_screen**
[2023-07-06 09:55:19] INFO: GTDB-Tk v2.3.0
[2023-07-06 09:55:19] INFO: gtdbtk classify --genome_dir genomes --align_dir align --out_dir classify -x gz --cpus 2 --skip_ani_screen
[2023-07-06 09:55:19] INFO: Using GTDB-Tk reference data version r214: /refdata/
[2023-07-06 09:55:20] TASK: Placing 2 archaeal genomes into reference tree with pplacer using 2 CPUs (be patient).
[2023-07-06 09:55:20] INFO: pplacer version: v1.1.alpha19-0-g807f6f3
[2023-07-06 09:59:33] INFO: Calculating RED values based on reference tree.
[2023-07-06 09:59:34] TASK: Traversing tree to determine classification method.
[2023-07-06 09:59:34] INFO: Completed 2 genomes in 0.00 seconds (15,448.63 genomes/second).
[2023-07-06 09:59:34] TASK: Calculating average nucleotide identity using FastANI (v1.32).
[2023-07-06 09:59:35] INFO: Completed 6 comparisons in 1.76 seconds (3.41 comparisons/second).
[2023-07-06 09:59:36] INFO: 2 genome(s) have been classified using FastANI and pplacer.
[2023-07-06 09:59:36] INFO: Note that Tk classification mode is insufficient for publication of new taxonomic designations. New designations should be based on one or more de novo trees, an example of which can be produced by Tk in de novo mode.
[2023-07-06 09:59:36] INFO: Done.
[user@cn3144 ~]$ **ls -lh classify**
[user@cn3144 ~]$ **cp -r /lscratch/$SLURM\_JOB\_ID/{identify,align,classify} .**
total 12K
drwxr-xr-x 3 user group 4.0K Jul  6 09:59 classify
lrwxrwxrwx 1 user group   32 Jul  6 09:59 gtdbtk.ar53.summary.tsv -> classify/gtdbtk.ar53.summary.tsv
-rw-r--r-- 1 user group 1.2K Jul  6 09:59 gtdbtk.json
-rw-r--r-- 1 user group 1.3K Jul  6 09:59 gtdbtk.log
-rw-r--r-- 1 user group    0 Jul  6 09:55 gtdbtk.warnings.log

```

A note from the documentation about the ANI screen:



> 
>  Starting with GTDB-Tk v2.2+, the classify\_wf and classify function
>  require an extra parameter to run: `--mash_db` or
>  `--skip_ani_screen`. With this new version of Tk, The first
>  stage of classify pipelines (classify\_wf and classify) is to compare all
>  user genomes to all reference genomes and annotate them, if possible, based
>  on ANI matches. Using the `--mash_db` option will indicate to
>  GTDB-Tk the path of the sketched Mash database require for ANI screening.
>  If no database are available ( i.e. this is the first time running classify
>  ), the `--mash_db` option will sketch a new Mash database that
>  can be used for subsequent calls. The `--skip_ani_screen` option
>  will skip the pre-screening step and classify all genomes similar to
>  previous versions of GTDB-Tk.1
> 


We provide a prebuilt mash sketch of the genomes in the reference data for
the classification step. This database was built with a kmer size of 16 and
sketch size of 5000 in `$GTDBTK_DB/mash/mash_16_5000.msh`. Since
this application is containerized this sketch is available as
`/refdata/mash/mash_16_5000.msh` in the container. Here is how the
last step would be run with the ANI screen using the pre-prepared mash database:



```

[user@cn3144 ~]$ **gtdbtk classify --genome\_dir genomes \
 --tmpdir /lscratch/$SLURM\_JOB\_ID/tmp \
 --align\_dir align \
 --out\_dir /lscratch/$SLURM\_JOB\_ID/classify \
 -x gz \
 --cpus $SLURM\_CPUS\_PER\_TASK \
 --mash\_db /refdata/mash/mash\_16\_5000.msh**
[2023-07-07 16:45:33] INFO: GTDB-Tk v2.3.2
[2023-07-07 16:45:33] INFO: gtdbtk classify --genome_dir genomes --align_dir align --out_dir classify -x gz --cpus 2 --mash_db /refdata/mash/mash_16_5000.msh
[2023-07-07 16:45:33] INFO: Using GTDB-Tk reference data version r214: /refdata/
[2023-07-07 16:45:34] INFO: Loading reference genomes.
[2023-07-07 16:45:34] INFO: Using Mash version 2.2.2
[2023-07-07 16:45:34] INFO: Creating Mash sketch file: classify/classify/ani_screen/intermediate_results/mash/gtdbtk.user_query_sketch.msh
[2023-07-07 16:45:34] INFO: Completed 2 genomes in 0.14 seconds (14.08 genomes/second).
[2023-07-07 16:45:34] INFO: Loading data from existing Mash sketch file: /refdata/mash/mash_16_5000.msh
[2023-07-07 16:45:41] INFO: Calculating Mash distances.
[2023-07-07 16:45:50] INFO: Calculating ANI with FastANI v1.32.
[2023-07-07 16:45:51] INFO: Completed 4 comparisons in 1.37 seconds (2.92 comparisons/second).
[2023-07-07 16:45:51] INFO: 2 genome(s) have been classified using the ANI pre-screening step.
[2023-07-07 16:45:51] INFO: Note that Tk classification mode is insufficient for publication of new taxonomic designations. New designations should be based on one or more de novo trees, an example of which can be produced by Tk in de novo mode.
[2023-07-07 16:45:51] INFO: Done.

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. gtdb-tk.sh). For example:



```

#!/bin/bash
set -e
module load gtdb-tk/2.3.0
mkdir /lscratch/$SLURM_JOB_ID/tmp
trap 'mv /lscratch/$SLURM_JOB_ID/classify classify_$SLURM_JOB_ID'
gtdbtk classify_wf --genome_dir genomes --out_dir /lscratch/$SLURM_JOB_ID/classify -x gz --skip_ani_screen --cpus $SLURM_CPUS_PER_TASK --tmpdir /lscratch/$SLURM_JOB_ID/tmp

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=2 --mem=50g --gres=lscratch:50 gtdb-tk.sh
```







