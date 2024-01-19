

document.querySelector('title').textContent = 'tetoolkit on Biowulf';
tetoolkit on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
 |


 A package for including transposable elements in differential enrichment
analysis of sequencing datasets.


### References:


* Ying Jin, O. H. Tam, E. Paniagua and M. Hammell. *TEtranscripts: a 
 package for including transposable elements in differential expression 
 analysis of RNA-seq datasets*. Bioinformatics 2015, 23:btv422.
 [Pubmed](http://www.ncbi.nlm.nih.gov/pubmed/26206304) | 
 PMC | 
 [Journal](http://bioinformatics.oxfordjournals.org/content/early/2015/08/06/bioinformatics.btv422.long)


Documentation
* [GitHub](https://github.com/mhammell-laboratory/tetoolkit)
* [Home page](http://hammelllab.labsites.cshl.edu/software/#TEToolkit)


Important Notes
* Module Name: tetoolkit (see [the modules page](/apps/modules.html) for more information)
* Example files in `$TETOOLKIT_TEST_DATA`
* Reference data in `/fdb/tetoolkit/`



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=10g --cpus-per-task=12 --gres=lscratch:50**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144]$ **module load tetoolkit**
[user@cn3144]$ **cp -rL $TETOOLKIT\_TEST\_DATA/01fastq .**
[user@cn3144]$ **cp $TETOOLKIT\_TEST\_DATA/GTF/dm6\_rmsk\_TE.gtf.gz .**
[user@cn3144]$ **gunzip dm6\_rmsk\_TE.gtf.gz**

```

This test data contains 4 untreated samples (GSM461176, GSM461177, GSM461178, GSM461182) and
3 samples with RNAi depleted *Pasilla* gene expression (GSM461179, GSM461180, GSM461181).
First, let's align it in a way that is compatible with TEtranscripts. This will take about 1h.



```

[user@cn3144]$ **module load STAR/2.7.3a**
[user@cn3144]$ **genome=/fdb/STAR\_indices/2.7.3a/UCSC/dm6/genes-75**
[user@cn3144]$ **mkdir -p 02aln**
[user@cn3144]$ **for fq in 01fastq/\*.fastq.gz ; do
 STAR --runThreadN $SLURM\_CPUS\_PER\_TASK \
 --genomeDir $genome \
 --readFilesCommand zcat \
 --outSAMtype BAM Unsorted \
 --winAnchorMultimapNmax 200 \
 --outFilterMultimapNmax 100 \
 --readFilesIn "$fq" \
 --outFileNamePrefix "./02aln/$(basename $fq .fastq.gz)"
 done**

```

Now we can run TEtoolkit on the alignments.



```

[user@cn3144]$ **untreated=(
 02aln/GSM461176Aligned.out.bam
 02aln/GSM461177Aligned.out.bam
 02aln/GSM461178Aligned.out.bam
 02aln/GSM461182Aligned.out.bam )**
[user@cn3144]$ **treated=(
 02aln/GSM461179Aligned.out.bam
 02aln/GSM461180Aligned.out.bam
 02aln/GSM461181Aligned.out.bam )**
[user@cn3144]$ **TEtranscripts --mode multi \
 --TE dm6\_rmsk\_TE.gtf \
 --GTF /fdb/STAR\_indices/2.7.3a/UCSC/dm6/genes.gtf \
 --minread 100 \
 --project test \
 -c "${untreated[@]}" \
 -t "${treated[@]}"**
...
[user@cn3144]$ **ls -1**
drwxr-xr-x 2 user group 4.0K Apr 10 08:50 01fastq
drwxr-xr-x 2 user group 4.0K Apr 10 09:38 02aln
-rw-r--r-- 1 user group 5.3M Apr 10 08:51 dm6_rmsk_TE.gtf
-rw-r--r-- 1 user group 521K Apr 10 14:13 test.cntTable
-rw-r--r-- 1 user group  735 Apr 10 14:13 test_DESeq2.R
-rw-r--r-- 1 user group 878K Apr 10 14:13 test_gene_TE_analysis.txt
-rw-r--r-- 1 user group  94K Apr 10 14:13 test_sigdiff_gene_TE.txt

[user@cn3144]$ **head -5 test\_sigdiff\_gene\_TE.txt**
baseMean        log2FoldChange  lfcSE   stat    pvalue  padj
AGBE    1668.50772310032        0.489686064951153       0.162938655205519       3.00534004244418        0.00265284088823233   0.0288848429674012
AGO2    8993.79284144723        -0.265578080875274      0.0892885493032728      -2.97438006270239       0.0029358119956592    0.0313958714926737
ATP7    1750.30129130752        0.526853871760357       0.161346748429612       3.26535165343106        0.00109328255464026   0.0139041285415587
Aats-cys        1256.3667600502 -0.670235513690024      0.105205551162389       -6.37072384759898       1.8813812974601e-10   1.13699601511619e-08

[user@cn3144]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. tetoolkit.sh), which uses the input file 'tetoolkit.in'. For example:



```

#! /bin/bash
# this is tetranscripts.sh

module load tetoolkit/2.1.4 || exit 1
untreated=(
     02aln/GSM461176Aligned.out.bam
     02aln/GSM461177Aligned.out.bam
     02aln/GSM461178Aligned.out.bam
     02aln/GSM461182Aligned.out.bam )
treated=(
     02aln/GSM461179Aligned.out.bam
     02aln/GSM461180Aligned.out.bam
     02aln/GSM461181Aligned.out.bam )
TEtranscripts --mode multi \
    --TE dm6_rmsk_TE.gtf \
    --GTF /fdb/STAR_indices/2.7.3a/UCSC/dm6/genes.gtf \
    --minread 100 \
    --project test \
    -c "${untreated[@]}" \
    -t "${treated[@]}"

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=2 --mem=10g tetoolkit.sh
```







