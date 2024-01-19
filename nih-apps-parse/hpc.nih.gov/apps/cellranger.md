

document.querySelector('title').textContent = 'cellranger on Biowulf';
cellranger on Biowulf


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


From the Cell Ranger manual:



> 
> Cell Ranger is a set of analysis pipelines that processes Chromium single cell
> 3’ RNA-seq output to align reads, generate gene-cell matrices and perform
> clustering and gene expression analysis. There are several pipelines:
> 
> * *cellranger mkfastq* wraps Illumina's bcl2fastq to correctly
>  demultiplex Chromium-prepared sequencing samples and to convert barcode
>  and read data to FASTQ files.
> * *cellranger count* takes FASTQ files from cellranger mkfastq and
>  performs alignment, filtering, and UMI counting. It uses the Chromium
>  cellular barcodes to generate gene-cell matrices and perform clustering
>  and gene expression analysis.
> * *cellranger aggr* aggregates results from cellranger count.
> * *cellranger reanalyze* takes feature-barcode matrices produced
>  by cellranger count or aggr and re-runs the dimensionality reduction, clustering,
>  and gene expression algorithms.
> * *cellranger multi* supports the anlysis of cell multiplexed data.
> 
> 
> Note that the command line interface has changed since version 1.1.
> 
> 
> 
> 
> These pipelines combine Chromium-specific algorithms with the widely used
> RNA-seq aligner STAR. Output is delivered in standard BAM, MEX, CSV, and HTML
> formats that are augmented with cellular information.
> 


Documentation
* [Home page + Manual](http://support.10xgenomics.com/single-cell/software/pipelines/latest/what-is-cell-ranger)


Important Notes
* Module Name: cellranger (see [the modules page](/apps/modules.html) for more information)
* cellranger can operate in local mode or 
 cluster mode. In both cases, the local part of the job will use
 multiple CPUs. Users have to specify the number of allocated CPUs and amount of memory
 with `--localcores=# --localmem=#` to cellranger.
* Please try local mode first and *only* use slurm mode if local mode does not produce results in a reasonable
 time frame. cellranger slurm mode tends to generate too many short jobs for moderate input sizes.
* cellranger may attempt to start more processes or open more files than the default limits
 on our compute nodes allow. If you encounter errors or strange results, you may have to raise these limits.
 See below for more deails.
* If using `cellranger multi` please allocate disk space on `lscratch`, assign the `TEMP` environment variable to `lscratch` and make sure the output is written to `lscratch`. Please email `staff@hpc.nih.gov` if you need help with this.
 * Reference data can be found in `$CELLRANGER_REF`. There are also environment variables for versioned subdirectories
 though their use is deprecated.
* Test data can be found in `$CELLRANGER_TEST_DATA`



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:


Copy the bcl format test data and run the demux pipeline



```

[user@biowulf]$ **sinteractive --cpus-per-task=6 --mem=35g**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load cellranger/5.0.0**
[user@cn3144 ~]$ **cp ${CELLRANGER\_TEST\_DATA:-none}/cellranger-tiny-bcl-1.2.0.tar.gz .**
[user@cn3144 ~]$ **cp ${CELLRANGER\_TEST\_DATA:-none}/cellranger-tiny-bcl-samplesheet-1.2.0.csv .**
[user@cn3144 ~]$ **tar -xzf cellranger-tiny-bcl-1.2.0.tar.gz**
[user@cn3144 ~]$ **cellranger mkfastq --run=cellranger-tiny-bcl-1.2.0 \
 --samplesheet=cellranger-tiny-bcl-samplesheet-1.2.0.csv \
 --localcores=$SLURM\_CPUS\_PER\_TASK \
 --localmem=34**
cellranger mkfastq (1.2.1)
Copyright (c) 2016 10x Genomics, Inc.  All rights reserved.
-------------------------------------------------------------------------------

Martian Runtime - 1.2.1 (2.1.2)
Running preflight checks (please wait)...
Checking run folder...
Checking RunInfo.xml...
Checking system environment...
Checking barcode whitelist...
Checking read specification...
Checking samplesheet specs...
2016-12-21 12:27:44 [runtime] (ready)           ID.H35KCBCXY.MAKE_FASTQS_CS.MAKE_FASTQS.PREPARE_SAMPLESHEET
[...snip...]
Outputs:
- Run QC metrics:        /spin1/users/user/test_data/cellranger/H35KCBCXY/outs/qc_summary.json
- FASTQ output folder:   /spin1/users/user/test_data/cellranger/H35KCBCXY/outs/fastq_path
- Interop output folder: /spin1/users/user/test_data/cellranger/H35KCBCXY/outs/interop_path
- Input samplesheet:     /spin1/users/user/test_data/cellranger/H35KCBCXY/outs/input_samplesheet.csv

Pipestance completed successfully!


```

Note that it is necessary to specify
`--localcores` and `--localmem`.


Cellranger may start an unreasonable number of processes or open too many
files. If you encounter errors that include



```

...
 self.pid = os.fork()
OSError: [Errno 11] Resource temporarily unavailable 

```

or see unexpected results despite specifying `--localcores` and
`--localmem`, you may have to raise the limit on the number of
processes and/or open files allowed in your batch script:



```

[user@cn3144 ~]$ **ulimit -u 10240 -n 16384**
[user@cn3144 ~]$ **cellranger mkfastq --run=cellranger-tiny-bcl-1.2.0 \
 --samplesheet=cellranger-tiny-bcl-samplesheet-1.2.0.csv \
 --localcores=$SLURM\_CPUS\_PER\_TASK \
 --localmem=34**

```

Generate counts per gene per cell



```

[user@cn3144 ~]$ **cellranger count --id s1 \
 --fastqs H35KCBCXY/outs/fastq\_path \
 --transcriptome=$CELLRANGER\_REF/refdata-gex-GRCh38-2020-A \
 --localcores=$SLURM\_CPUS\_PER\_TASK \
 --chemistry=SC3Pv2 \
 --localmem=34**
cellranger count (1.2.1)
Copyright (c) 2016 10x Genomics, Inc.  All rights reserved.
-------------------------------------------------------------------------------

Martian Runtime - 1.2.1 (2.1.2)
Running preflight checks (please wait)...
Checking sample info...
Checking FASTQ folder...
[...snip...]
Outputs:
- Run summary HTML:                      /spin1/users/user/test_data/cellranger/s1/outs/web_summary.html
- Run summary CSV:                       /spin1/users/user/test_data/cellranger/s1/outs/metrics_summary.csv
- BAM:                                   /spin1/users/user/test_data/cellranger/s1/outs/possorted_genome_bam.bam
- BAM index:                             /spin1/users/user/test_data/cellranger/s1/outs/possorted_genome_bam.bam.bai
- Filtered gene-barcode matrices MEX:    /spin1/users/user/test_data/cellranger/s1/outs/filtered_gene_bc_matrices
- Filtered gene-barcode matrices HDF5:   /spin1/users/user/test_data/cellranger/s1/outs/filtered_gene_bc_matrices_h5.h5
- Unfiltered gene-barcode matrices MEX:  /spin1/users/user/test_data/cellranger/s1/outs/raw_gene_bc_matrices
- Unfiltered gene-barcode matrices HDF5: /spin1/users/user/test_data/cellranger/s1/outs/raw_gene_bc_matrices_h5.h5
- Secondary analysis output CSV:         /spin1/users/user/test_data/cellranger/s1/outs/analysis
- Per-molecule read information:         /spin1/users/user/test_data/cellranger/s1/outs/molecule_info.h5

Pipestance completed successfully!

Saving pipestance info to s1/s1.mri.tgz

```

The same job could also be run in cluster mode where pipeline tasks
are submitted as batch jobs. This can be done by setting jobmode to slurm
and limiting the max. number of concurrent jobs:



```

[user@cn3144 ~]$ **cellranger count --id s1 \
 --fastqs H35KCBCXY/outs/fastq\_path \
 --transcriptome=$CELLRANGER\_REF/refdata-gex-GRCh38-2020-A \
 --chemistry=SC3Pv2 \
 --localcores=$SLURM\_CPUS\_PER\_TASK \
 --localmem=34 \
 --jobmode=slurm --maxjobs=10**

```

If running in slurm mode, it may be necessary to add `--jobinterval=3000` if encountering
errors mentioning empty batch scripts.


Though in the case of this small example this actually results in
a longer overall runtime. Even when running in cluster mode, please run
the main pipeline in an sinteractive session or as a batch job itself.


Don't forget to close the interactive session when done



```

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```

### Interesting bits


As per [10x genomics](https://www.10xgenomics.com/resources/analysis-guides/tutorial-navigating-10x-barcoded-bam-files#:~:text=To%20use%20this%20samtools%20you%20can%20run%20the,file%20and%20the%20time%20required%20to%20sort%20it.) the `xf` auxillary field in the
output bam file can be used to filter a bam file to include only the representative of a UMI confidently mapped to a feature and/or
transcriptome. This can be done with versions of samtools that include the `-e` flag to filter by an expression:



```

[user@cn3144 ~]$ **module load samtools/1.17**
[user@cn3144 ~]$ **samtools view -@6 -e '[xf]==25' -o filtered.bam s1/outs/possorted\_genome\_bam.bam**

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. cellranger.sh), which uses the input file 'cellranger.in'. For example:



```

#! /bin/bash
module load cellranger/5.0.0 || exit 1
## uncomment the following line if encountering 'resource unavailable' errors
## despite using --localcores and --localmem
# ulimit -u 4096
cellranger mkfastq --run=llranger-tiny-bcl-1.2.0 \
        --samplesheet=cellranger-tiny-bcl-samplesheet-1.2.0.csv \
        --localcores=$SLURM_CPUS_PER_TASK \
        --localmem=34
cellranger count --id s1 \
        --fastqs H35KCBCXY/outs/fastq_path \
        --transcriptome=$CELLRANGER_REF/refdata-gex-GRCh38-2020-A \
        --localcores=$SLURM_CPUS_PER_TASK \
        --localmem=34 \
        --chemistry=SC3Pv2 \
        --jobmode=slurm --maxjobs=20

```

Again, please remember to include `--localcoes` and `--localmem`.


Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=12 --mem=35g cellranger.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. cellranger.swarm). For example:



```

cellranger mkfastq --run=./run1 --localcores=$SLURM_CPUS_PER_TASK --localmem=34
cellranger mkfastq --run=./run2 --localcores=$SLURM_CPUS_PER_TASK --localmem=34
cellranger mkfastq --run=./run3 --localcores=$SLURM_CPUS_PER_TASK --localmem=34

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f cellranger.swarm -g 35 -t 12 --module cellranger
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g #  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t #  Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module cellranger  Loads the cellranger module for each subjob in the swarm 
 | |
 | |
 | |








