

document.querySelector('title').textContent = 'spaceranger on Biowulf';
spaceranger on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
 |




 From the 10x spaceranger documentation:
 



>  Space Ranger is a set of analysis pipelines that process
>  Visium spatial RNA-seq output and brightfield microscope images in
>  order to detect tissue, align reads, generate feature-spot matrices,
>  perform clustering and gene expression analysis, and place spots in
>  spatial context on the slide image 




Documentation
* Spaceranger[Manual](https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/what-is-space-ranger)


Important Notes
* Module Name: spaceranger (see [the modules page](/apps/modules.html) 
 for more information)
* spaceranger can operate in local mode or 
 cluster mode. In both cases, the local part of the job will use
 multiple CPUs. Users have to specify the number of allocated CPUs and amount of memory
 with `--localcores=# --localmem=#` to spaceranger.
* spaceranger may attempt to start more processes or open more files than the default limits
 on our compute nodes allow. If you encounter errors or strange results, you may have to raise these limits.
 See below for more deails.
* Test data can be found in `$SPACERANGER_TEST_DATA`
* Reference data can be found in `$SPACERANGER_REF`



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.


First step is to generate fastq files from bcl using the spaceranger
adapter. This is similar to all the 10x tools. For this example we will
use the tiny bcl data set used in the 10x genomics manual.



```

[user@biowulf]$ **sinteractive --cpus-per-task=16 --mem=38g**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144]$ **ml spaceranger**
[user@cn3144]$ **cp -r $SPACERANGER\_TEST\_DATA/mkfastq data1**
[user@cn3144]$ **ls -lh data1**
total 1.3G
-rw-r--r-- 1 user group 1.3G Dec 11 15:53 spaceranger-tiny-bcl-1.0.0.tar.gz
-rw-r--r-- 1 user group  552 Dec 11 15:53 spaceranger-tiny-bcl-samplesheet-1.0.0.csv
-rw-r--r-- 1 user group   41 Dec 11 15:53 spaceranger-tiny-bcl-simple-1.0.0.csv
[user@cn3144]$ **cat data1/spaceranger-tiny-bcl-simple-1.0.0.csv**
Lane,Sample,Index
1,test_sample,SI-TT-D9
[user@cn3144]$ **cd data1 && tar -xzf spaceranger-tiny-bcl-1.0.0.tar.gz && cd ..**
[user@cn3144]$ **spaceranger mkfastq --id=tiny-bcl \
 --run=data1/spaceranger-tiny-bcl-1.0.0 \
 --csv=data1/spaceranger-tiny-bcl-simple-1.0.0.csv**
spaceranger mkfastq (spaceranger-1.2.2)
Copyright (c) 2021 10x Genomics, Inc.  All rights reserved.
-------------------------------------------------------------------------------

The --qc option is deprecated and will be removed in a future version.
Most qc metrics can be found in the web summary.
Martian Runtime - v4.0.2

Running preflight checks (please wait)...
Checking run folder...
Checking RunInfo.xml...
Checking system environment...
Emitting run information...
Checking read specification...
Checking samplesheet specs...
[...snip...]
Pipestance completed successfully!


```

Next - counting reads for the capture areas. Note that for this step
we will use already demultiplexted data from a mouse brain section obtained
from [10X genomics](https://support.10xgenomics.com/spatial-gene-expression/datasets/1.0.0/V1_Mouse_Brain_Sagittal_Posterior).



```

[user@cn3144]$ **cp -rL $SPACERANGER\_TEST\_DATA/count data2** 
[user@cn3144]$  **spaceranger count --id=test \
 --transcriptome=${SPACERANGER\_REF}/refdata-gex-mm10-2020-A \
 --fastqs=data2/V1\_Mouse\_Brain\_Sagittal\_Posterior\_Section\_1\_fastqs \
 --sample=V1\_Mouse\_Brain\_Sagittal\_Posterior\_Section\_1 \
 --image=data2/V1\_Mouse\_Brain\_Sagittal\_Posterior\_image.tif \
 --slide=V19L29-035 \
 --area=A1 --localcores=$SLURM\_CPUS\_PER\_TASK --localmem=37**
...many hours later...
Waiting 6 seconds for UI to do final refresh.
Pipestance completed successfully!

2019-12-12 21:27:20 Shutting down.
Saving pipestance info to "test/test.mri.tgz"

```

Note that it is necessary to specify
`--localcores` and `--localmem`.


Spaceranger may start an unreasonable number of processes or open too many
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

The same job could also be run in cluster mode where pipeline tasks
are submitted as batch jobs. This can be done by setting jobmode to slurm
and limiting the max. number of concurrent jobs:



```

[user@cn3144]$  **spaceranger count --id=test \
 --transcriptome=${SPACERANGER\_REF}/refdata-gex-mm10-2020-A \
 --fastqs=data2/V1\_Mouse\_Brain\_Sagittal\_Posterior\_Section\_1\_fastqs \
 --sample=V1\_Mouse\_Brain\_Sagittal\_Posterior\_Section\_1 \
 --image=data2/V1\_Mouse\_Brain\_Sagittal\_Posterior\_image.tif \
 --slide=V19L29-035 \
 --area=A1 --localcores=$SLURM\_CPUS\_PER\_TASK --localmem=37 \
 --jobmode=slurm --maxjobs=20**
...many hours later...
Waiting 6 seconds for UI to do final refresh.
Pipestance completed successfully!

2019-12-12 21:27:20 Shutting down.
Saving pipestance info to "test/test.mri.tgz"
[user@cn3144]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. spaceranger.sh), which uses the input file 'spaceranger.in'. For example:



```

#!/bin/bash
module load spaceranger
spaceranger count --id=test \
    --transcriptome=${SPACERANGER_REF}/refdata-gex-mm10-2020-A \
    --fastqs=data2/V1_Mouse_Brain_Sagittal_Posterior_Section_1_fastqs \
    --sample=V1_Mouse_Brain_Sagittal_Posterior_Section_1 \
    --image=data2/V1_Mouse_Brain_Sagittal_Posterior_image.tif \
    --slide=V19L29-035 \
    --area=A1 --localcores=$SLURM_CPUS_PER_TASK --localmem=37 \
    --jobmode=slurm --maxjobs=20

```

Again, please remember to include `--localcores` and `--localmem`


Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=8 --mem=38g spaceranger.sh
```







