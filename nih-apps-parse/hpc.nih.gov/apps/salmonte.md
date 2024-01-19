

document.querySelector('title').textContent = 'salmonte on Biowulf';
salmonte on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
 |


 SalmonTE is an ultra-Fast and Scalable Quantification Pipeline of Transpose
Element (TE) Abundances. It is based on snakemake, salmon and R. Note that SalmonTE
packages its own version of salmon.


### References:


* H. H. Jeong, H. K. Yalamanchili, C. Guo C, J. M. Shulman, and Z. Liu
*An ultra-fast and scalable quantification pipeline for transposable elements from next generation sequencing data.* Pac Symp Biocomput. 23:168-179 (2018).
[PubMed](https://www.ncbi.nlm.nih.gov/pubmed/29218879) | 
[Journal](https://www.worldscientific.com/doi/abs/10.1142/9789813235533_0016)


Documentation
* SalmonTE on [GitHub](https://github.com/LiuzLab/SalmonTE)


Important Notes
* Module Name: salmonte (see [the modules page](/apps/modules.html)
for more information)
* `SalmonTE.py quant` is multithreaded. Please match the number of
threads to the number of allocated CPUs
* Example files in `$SALMONTE_TEST_DATA`
* Reference data in /fdb/salmonte/



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=6g --cpus-per-task=4 --gres=lscratch:10**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144]$ **cd /lscratch/$SLURM\_JOB\_ID**
[user@cn3144]$ **module load salmonte**
[user@cn3144]$ **cp -r ${SALMONTE\_TEST\_DATA:-none}/data .**
[user@cn3144]$ **ls -lh data**
total 5.0M
-rw-rw-r-- 1 user group 634K Nov 11 10:13 CTRL_1_R1.fastq
-rw-rw-r-- 1 user group 634K Nov 11 10:13 CTRL_1_R2.fastq
-rw-rw-r-- 1 user group 634K Nov 11 10:13 CTRL_2_R1.fastq
-rw-rw-r-- 1 user group 634K Nov 11 10:13 CTRL_2_R2.fastq
-rw-rw-r-- 1 user group 634K Nov 11 10:13 TARDBP_1_R1.fastq
-rw-rw-r-- 1 user group 634K Nov 11 10:13 TARDBP_1_R2.fastq
-rw-rw-r-- 1 user group 634K Nov 11 10:13 TARDBP_2_R1.fastq
-rw-rw-r-- 1 user group 634K Nov 11 10:13 TARDBP_2_R2.fastq
[user@cn3144]$ **SalmonTE.py quant --reference=hs --outpath=quant\_out \
 --num\_threads=$SLURM\_CPUS\_PER\_TASK --exprtype=count data**
2019-11-11 10:19:30,550 Starting quantification mode
2019-11-11 10:19:30,550 Collecting FASTQ files...
2019-11-11 10:19:30,553 The input dataset is considered as a paired-ends dataset.
2019-11-11 10:19:30,553 Collected 4 FASTQ files.
2019-11-11 10:19:30,553 Quantification has been finished.
2019-11-11 10:19:30,553 Running Salmon using Snakemake
...
[user@cn3144]$ **ls -lh quant\_out**
total 68K
-rw-rw-r-- 1 user group  23K Nov 11 10:19 clades.csv
-rw-rw-r-- 1 user group   63 Nov 11 10:19 condition.csv
drwxrwxr-x 5 user group 4.0K Nov 11 10:19 CTRL_1
drwxrwxr-x 5 user group 4.0K Nov 11 10:19 CTRL_2
-rw-rw-r-- 1 user group  17K Nov 11 10:19 EXPR.csv
-rw-rw-r-- 1 user group  161 Nov 11 10:19 MAPPING_INFO.csv
drwxrwxr-x 5 user group 4.0K Nov 11 10:19 TARDBP_1
drwxrwxr-x 5 user group 4.0K Nov 11 10:19 TARDBP_2

```

Notes:


1. SalmonTE tries to determine if data is paired end by comparing sequence
 ids. This can fail with certain id lines like, for example, those
 produced by the sra-toolkit. In other words, it's a fragile
 autodetection and fastq files may have to be reformatted to make it
 work correctly.
2. SalmonTE supports compressed fastq files


Before running the differential expression test, it is necessary
to update the file quant\_out/condition.csv to include your experimental
conditions.



```

[user@cn3144]$ **mv quant\_out/condition.csv quant\_out/condition.csv.orig**
[user@cn3144]$ **cat <<EOF > quant\_out/condition.csv**
SampleID,condition
TARDBP_1,treatment
CTRL_1,control
TARDBP_2,treatment
CTRL_2,control
**EOF**
[user@cn3144]$ ### or just edit the condition.csv file with your favorite text editor
[user@cn3144]$ **SalmonTE.py test --inpath=quant\_out --outpath=test\_out \
 --tabletype=csv --figtype=png --analysis\_type=DE \
 --conditions=control,treatment**
[user@cn3144]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. salmonte.sh), which uses the input file 'salmonte.in'. For example:



```

#!/bin/bash
module load salmonte/0.4 || exit 1
SalmonTE.py quant  --reference=hs --outpath=all_quant_out \
    --num_threads=$SLURM_CPUS_PER_TASK --exprtype=count data

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=8 --mem=10g salmonte.sh
```







