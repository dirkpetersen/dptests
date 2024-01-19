

document.querySelector('title').textContent = 'Cogent NGS Analysis Pipeline on Biowulf';

 .hl { background-color: #ffff99; }

Cogent NGS Analysis Pipeline on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
 |



From the TaKaRa manual page:



> 
> 
>  Cogent NGS Analysis Pipeline (CogentAP) is bioinformatic software for analyzing
>  RNA-seq NGS data generated using the following systems or kits:
> 
>  * ICELL8 cx Single-Cell System or the ICELL8 Single-Cell System on the single-cell full-length transcriptome (SMART-Seq ICELL8 workflow)
> * ICELL8 cx Single-Cell System or the ICELL8 Single-Cell System on the single-cell differential expression (3′ DE or 5′ DE) workflows (ICELL8 3′ DE or ICELL8 TCR)
> * SMARTer Stranded Total RNA-Seq Kit v3 - Pico Input Mammalian
> 
> 
> 
>  The program takes input data from sequencing and outputs an HTML
>  report, with results typical to single-cell analysis, plus other
>  files, such as a gene matrix, to continue further analysis. R data
>  object with pre-computed results based on recommended parameters
>  are also output. Either the standard output files or the R data
>  object can serve as input for Cogent NGS Discovery Software
>  (CogentDS), another bioinformatic software package provided by
>  Takara Bio.
> 
>  CogentAP software is written in Python and can be run either via a
>  GUI or command-line interface.
> 
>  



Documentation
* CogentAP [manual](https://www.takarabio.com/learning-centers/next-generation-sequencing/bioinformatics-resources/cogent-ngs-analysis-pipeline/cogent-ngs-analysis-pipeline-v10-user-manual)


Important Notes
* Module Name: cogentap (see [the modules page](/apps/modules.html) for more information)
* Multithreaded/singlethreaded/MPI...
* Example files in `$COGENTAP_TEST_DATA`
* Reference data is stored in `/fdb/cogentap` and linked into the
 install directory at the expected path



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --cpus-per-task=8 --mem=45g --gres=lscratch:20**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144]$ **cd /lscratch/$SLURM\_JOB\_ID**
[user@cn3144]$ **module load cogentap**
[user@cn3144]$ **cogent --help**
usage: cogent

Script to perform NGS analysis. Please see helps of each command for details.

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         Show version number

commands:
  {add_genome,demux,analyze}
    add_genome          Build a genome with preferred STAR parameters.
    demux               De-multiplex barcoded reads from sequence data stored in FASTQ files.
    analyze             Perform counting analysis for exons and genes by fastq input data.

[user@cn3144]$ **cp -r ${COGENTAP\_TEST\_DATA} .**
[user@cn3144]$ **cogent demux \
 -i test/test\_FL\_R1.fastq.gz \
 -p test/test\_FL\_R2.fastq.gz \
 --barcodes\_file test/99999\_CogentAP\_test\_selected\_WellList.TXT \
 -t ICELL8\_FLA \
 -o out \
 -n $SLURM\_CPUS\_PER\_TASK**
###
### cogent 1.0
###
[user@cn3144]$ **cogent analyze \
 -i out/out\_demuxed\_R1.fastq \
 -p out/out\_demuxed\_R2.fastq \
 -g hg38 \
 -o out/analysis \
 -n $SLURM\_CPUS\_PER\_TASK \
 -d out/out\_counts\_all.csv \
 -t ICELL8\_FLA**

###
### cogent ≥1.5.0 - see the official manual for more differences
###
[user@cn3144]$ **cogent analyze \
 -i out \
 -g hg38 \
 -o out/analysis \
 --threads $SLURM\_CPUS\_PER\_TASK \
 -t ICELL8\_FLA**

[user@cn3144]$ **tree out**
out
├ [user   4.0K]  analysis
│   ├ [user   5.2K]  analysis_analyzer.log
│   ├ [user   2.1M]  analysis_genematrix.csv
│   ├ [user   1.0K]  analysis_stats.csv
│   ├ [user   4.0K]  cogent_ds
│   │   ├ [user   1.8M]  CogentDS.analysis.rda
│   │   ├ [user   214K]  CogentDS.boxplot.png
│   │   ├ [user   3.4K]  CogentDS.cogent_ds.log
│   │   ├ [user     70]  CogentDS.cor_stats.csv
│   │   ├ [user   164K]  CogentDS.heatmap.png
│   │   ├ [user   1.7M]  CogentDS.report.html
│   │   └ [user   157K]  CogentDS.UMAP.png
│   ├ [user   4.0K]  extras
│   │   ├ [user   2.1M]  analysis_incl_introns_genematrix.csv
│   │   ├ [user   1.0K]  analysis_incl_introns_stats.csv
│   │   └ [user   3.9M]  gene_info_incl_introns.csv
│   ├ [user   3.8M]  gene_info.csv
│   └ [user   4.0K]  work
│       ├ [user    14M]  analysis.Aligned.out.bam
│       ├ [user   2.0K]  analysis.Log.final.out
│       ├ [user   488K]  analysis.SJ.out.tab
│       └ [user     39]  mito.gtf
├ [user    257]  out_counts_all.csv
├ [user    20M]  out_demuxed_R1.fastq
├ [user    20M]  out_demuxed_R2.fastq
└ [user   1.3K]  out_demuxer.log

[user@cn3144]$ **mv out /data/$USER/**
[user@cn3144]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. cogentap.sh) similar to the following:



```

#!/bin/bash
module load cogentap/1.0
cd /lscratch/$SLURM_JOB_ID || exit 1
module load cogentap
cp -r ${COGENTAP_TEST_DATA:-none} .
cogent demux \
        -i test/test_FL_R1.fastq.gz \
        -p test/test_FL_R2.fastq.gz \
        --barcodes_file test/99999_CogentAP_test_selected_WellList.TXT \
        -t ICELL8_FLA \
        -o out \
        -n $SLURM_CPUS_PER_TASK
cogent analyze \
        -i out/out_demuxed_R1.fastq \
        -p out/out_demuxed_R2.fastq \
        -g hg38 \
        -o out/analysis \
        -n $SLURM_CPUS_PER_TASK \
        -d out/out_counts_all.csv \
        -t ICELL8_FLA

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=8 --mem=30g cogentap.sh
```







