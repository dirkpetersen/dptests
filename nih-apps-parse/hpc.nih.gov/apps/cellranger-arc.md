

document.querySelector('title').textContent = 'cellranger-arc on Biowulf';
cellranger-arc on Biowulf


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


From the Cell Ranger Arc manual:



>  
>  Cell Ranger ARC is a set of analysis pipelines that process
>  Chromium Single Cell Multiome ATAC + Gene Expression sequencing data to
>  generate a variety of analyses pertaining to gene expression, chromatin
>  accessibility and their linkage. Furthermore, since the ATAC and gene
>  expression measurements are on the very same cell, we are able to perform
>  analyses that link chromatin accessibility and gene expression.
> 
> * *cellranger-arc mkfastq* demultiplexes raw base call (BCL)
>  files generated by Illumina sequencers into FASTQ files. It is a wrapper
>  around Illumina's bcl2fastq, with additional useful features that are
>  specific to 10x libraries and a simplified sample sheet format. The same
>  command can be used to demultiplex both ATAC and GEX flow cells.
> * *cellranger-arc count* takes FASTQ files from cellranger-arc
>  mkfastq and performs alignment, filtering, barcode counting, peak calling
>  and counting of both ATAC and GEX molecules. Furthermore, it uses the
>  Chromium cellular barcodes to generate feature-barcode matrices, perform
>  dimensionality reduction, determine clusters, perform differential analysis
>  on clusters and identify linkages between peaks and genes. The count
>  pipeline can take input from multiple sequencing runs on the same GEM
>  well.
> 
> 
> 
> 
> These pipelines combine Chromium-specific algorithms with the widely used
> RNA-seq aligner STAR. Output is delivered in standard BAM, MEX, CSV, HDF5, and HTML
> formats that are augmented with cellular information and a .cloupe file for use
> with the Loupe browser.
> 


Documentation
* [Home page + Manual](https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/pipelines/latest/what-is-cell-ranger-arc)


Important Notes
* Module Name: cellranger-arc (see [the modules page](/apps/modules.html) for more information)
* cellranger can operate in local mode or 
 cluster mode. In both cases, the local part of the job will use
 multiple CPUs. Users have to specify the number of allocated CPUs and amount of memory
 with `--localcores=# --localmem=#` to cellranger.
* cellranger may attempt to start more processes or open more files than the default limits
 on our compute nodes allow. If you encounter errors or strange results, you may have to raise these limits.
 See below for more deails.
* Reference data can be found in `$CELLRANGER_ARC_REF`
* Test data can be found in `$CELLRANGER_ARC_TEST_DATA`



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

[user@cn3144 ~]$ **module load cellranger-arc**
[user@cn3144 ~]$ **cp ${CELLRANGER\_ARC\_TEST\_DATA:-none}/\* .**
[user@cn3144 ~]$ **tar -xzf cellranger-arc-tiny-bcl-atac-1.0.0.tar.gz**
[user@cn3144 ~]$ **tar -xzf cellranger-arc-tiny-bcl-gex-1.0.0.tar.gz**
[user@cn3144 ~]$ # demultiplex the ATAC flowcell
[user@cn3144 ~]$ **cellranger-arc mkfastq --id=tiny-bcl-atac \
 --csv=cellranger-arc-tiny-bcl-atac-simple-1.0.0.csv \
 --run=cellranger-arc-tiny-bcl-atac-1.0.0 \
 --localcores=$SLURM\_CPUS\_PER\_TASK \
 --localmem=34**
cellranger-arc mkfastq (cellranger-arc-1.0.0)
Copyright (c) 2020 10x Genomics, Inc.  All rights reserved.
-------------------------------------------------------------------------------

Martian Runtime - v4.0.1
Serving UI at http://cn1038:38947?auth=OKNOP2vgDnOXFmUoe1dK2y4rCKFuNCkYs16KtaTMqfw

Running preflight checks (please wait)...
Checking run folder...
Checking RunInfo.xml...
...
Pipestance completed successfully!

2020-10-06 20:30:54 Shutting down.
Saving pipestance info to "tiny-bcl-atac/tiny-bcl-atac.mri.tgz"

[user@cn3144 ~]$ # demultiplex the GEX flowcell
[user@cn3144 ~]$ **cellranger-arc mkfastq --id=tiny-bcl-gex \
 --csv=cellranger-arc-tiny-bcl-gex-simple-1.0.0.csv \
 --run=cellranger-arc-tiny-bcl-gex-1.0.0 \
 --localcores=$SLURM\_CPUS\_PER\_TASK \
 --localmem=34**
...


```

Note that it is necessary to specify
`--localcores` and `--localmem`.


Cellranger Arc may start an unreasonable number of processes or open too many
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

```

If running in slurm mode, it may be necessary to add `--jobinterval=3000` if encountering
errors mentioning empty batch scripts.


Generate counts per gene per cell



```

[user@cn3144 ~]$ **cat > gex\_atac.csv <<\_\_EOF\_\_**
fastqs,sample,library_type
tiny-bcl-gex/outs/fastq_path/,test_sample_gex,Gene Expression
tiny-bcl-atac/outs/fastq_path/,test_sample_atac,Chromatin Accessibility
**\_\_EOF\_\_**
[user@cn3144 ~]$ **cellranger-arc count \
 --id=sample123 \
 --reference=$CELLRANGER\_ARC\_REF/refdata-cellranger-arc-GRCh38-2020-A \
 --libraries=gex\_atac.csv \
 --localcores=$SLURM\_CPUS\_PER\_TASK \
 --localmem=32**
### note: example data currently fails at this step. waiting for reply from 10x support

```

The same job could also be run in cluster mode where pipeline tasks
are submitted as batch jobs. This can be done by setting jobmode to slurm
and limiting the max. number of concurrent jobs:



```

[user@cn3144 ~]$ **cellranger-arc count \
 --jobmode=slurm --maxjobs=10 \
 --id=sample123 \
 --reference=$CELLRANGER\_ARC\_REF/refdata-cellranger-arc-GRCh38-2020-A \
 --libraries=gex\_atac.csv \
 --localcores=$SLURM\_CPUS\_PER\_TASK \
 --localmem=32**

```

Don't forget to close the interactive session when done



```

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. cellranger.sh), which uses the input file 'cellranger.in'. For example:



```

#! /bin/bash
module load cellranger || exit 1
## uncomment the following line if encountering 'resource unavailable' errors
## despite using --localcores and --localmem
# ulimit -u 4096
cellranger-arc mkfastq --id=tiny-bcl-atac \
     --csv=cellranger-arc-tiny-bcl-atac-simple-1.0.0.csv \
     --run=cellranger-arc-tiny-bcl-atac-1.0.0 \
     --localcores=$SLURM_CPUS_PER_TASK \
     --localmem=34
cellranger-arc mkfastq --id=tiny-bcl-gex \
     --csv=cellranger-arc-tiny-bcl-gex-simple-1.0.0.csv \
     --run=cellranger-arc-tiny-bcl-gex-1.0.0 \
     --localcores=$SLURM_CPUS_PER_TASK \
     --localmem=34

cat > gex_atac.csv <<__EOF__
fastqs,sample,library_type
tiny-bcl-gex/outs/fastq_path/,test_sample_gex,Gene Expression
tiny-bcl-atac/outs/fastq_path/,test_sample_atac,Chromatin Accessibility
__EOF__
cellranger-arc count \
      --id=sample123 \
      --reference=$CELLRANGER_ARC_REF/refdata-cellranger-arc-GRCh38-2020-A \
      --libraries=gex_atac.csv \
      --localcores=6 \
      --localmem=32

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

cellranger-arc mkfastq --run=./run1 --localcores=$SLURM_CPUS_PER_TASK --localmem=34
cellranger-arc mkfastq --run=./run2 --localcores=$SLURM_CPUS_PER_TASK --localmem=34
cellranger-arc mkfastq --run=./run3 --localcores=$SLURM_CPUS_PER_TASK --localmem=34

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








