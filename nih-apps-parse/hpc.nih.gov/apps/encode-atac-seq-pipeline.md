

document.querySelector('title').textContent = 'encode-atac-seq-pipeline on Biowulf';
encode-atac-seq-pipeline on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
 |




From the Encode documentation:




>  
>  This pipeline is designed for automated end-to-end quality control
>  and processing of ATAC-seq or DNase-seq data. The pipeline can be run on
>  compute clusters with job submission engines or stand alone machines. It
>  inherently makes uses of parallelized/distributed computing. 
> 



Documentation
* [Description](https://www.encodeproject.org/atac-seq) on the Encode site
* On [GitHub](https://github.com/ENCODE-DCC/atac-seq-pipeline)
* [Input Json](https://github.com/ENCODE-DCC/atac-seq-pipeline/blob/master/docs/input.md)
 description
* [Output](https://github.com/ENCODE-DCC/atac-seq-pipeline/blob/master/docs/output.md)
 description


Important Notes
* Module Name: encode-atac-seq-pipeline (see [the modules page](/apps/modules.html) for more information)
* This pipeline has undergone frequent updates with backwards incompatible
 changes to execution and configuration. If in doubt please refer back to the
 upstream documentation.
* For local runs the CPU and memory consumption varies over time and its magnitude
 depends on the number of replicates and size of input data. For the example below (2 fastq files per replicaate,
 2 replicates), 10-12 CPUs and 16GB of memory were sufficient.
* As of version 2.0, we no longer provide a Slurm backend configuration. We recommend use of this pipeline in 'local'
 backend mode. If you believe you have a use for Slurm backend execution, please get in touch with HPC staff.
* Environment variables set (not all of them are set in each version due to significant
changes in the pipeline):
	+ `$EASP_BACKEND_CONF`: configuration for local backend
	+ `$EASP_WFOPTS`: singularity backend opts (Versions > 1.0 only)
	+ `$EASP_WDL`: WDL file defining the workflow
	+ `$EASP_TEST_DATA`: Input data for example below* Reference data in /fdb/encode-atac-seq-pipeline/<version>



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:


A note about resource allocation:


* Set the number of concurrent tasks (NUM\_CONCURRENT\_TASKS) to the number of replicates
* Set `"atac.bowtie2_cpu"` in the input json to the number of CPUs you want
 bowtie2 to use. Usually 8 or so.
* Allocate `NUM_CONCURRENT_TASKS * atac.bowtie2_cpu` CPUs
* Allocate `20GB * NUM_CONCURRENT_TASKS` for big samples and 
 `10GB * NUM_CONCURRENT_TASKS` for small samples



WDL based workflows need a json file to define input and settings for a workflow run. In
this example, we will use the 50nt data from ENCODE sample 
[ENCSR889WQX](https://www.encodeproject.org/experiments/ENCSR889WQX/) (mouse frontal
cortex). This includes 2 fastq files for each of 2 replicates.




* 2.1.0





```


[user@biowulf]$ **sinteractive --cpus-per-task=8 --mem=20g --gres=lscratch:30**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144]$ **wd=$PWD**  # so we can copy results back later
[user@cn3144]$ **cd /lscratch/$SLURM\_JOB\_ID**
[user@cn3144]$ **module load encode-atac-seq-pipeline/1.0**
[user@cn3144]$ **cp -Lr ${EASP\_TEST\_DATA:-none}/\* .**
[user@cn3144]$ **tree**
.
|-- [user    738]  ENCSR889WQX.json.1.0
|-- [user    870]  ENCSR889WQX.json.1.1.7
`-- [user   4.0K]  input
    |-- [user   4.0K]  rep1
    |   |-- [user   213M]  ENCFF439VSY_5M.fastq.gz
    |   `-- [user   213M]  ENCFF683IQS_5M.fastq.gz
    `-- [user   4.0K]  rep2
        |-- [user   210M]  ENCFF463QCX_5M.fastq.gz
        `-- [user   211M]  ENCFF992TSA_5M.fastq.gz

[user@cn3144]$ **cat ENCSR889WQX.json.1.0**
{
    "atac.pipeline_type" : "atac",
    "atac.genome_tsv" : "/fdb/encode-atac-seq-pipeline/v1/mm10/mm10.tsv",
    "atac.fastqs" : [
        [
            ["input/rep1/ENCFF683IQS_5M.fastq.gz"],
            ["input/rep1/ENCFF439VSY_5M.fastq.gz"]
        ],
        [
            ["input/rep2/ENCFF992TSA_5M.fastq.gz"],
            ["input/rep2/ENCFF463QCX_5M.fastq.gz"]
        ]
    ],
    "atac.paired_end" : false,
    "atac.multimapping" : 4,
    "atac.trim_adapter.auto_detect_adapter" : true,
    "atac.smooth_win" : 73,
    "atac.enable_idr" : true,
    "atac.idr_thresh" : 0.05,
    "atac.qc_report.name" : "ENCSR889WQX (subsampled to 5M reads)",
    "atac.qc_report.desc" : "ATAC-seq on Mus musculus C57BL/6 frontal cortex adult"
}
        
```

In this example the pipeline will only be run locally - i.e. it will not submit
 tasks as slurm jobs



```

[user@cn3144]$ **java -Dconfig.file=$EASP\_BACKEND\_CONF \
 -Dbackend.default=Local \
 -jar $CROMWELL\_JAR run -i ENCSR889WQX.json.1.0 $EASP\_WDL**
[...much output...]
[user@cn3144]$ **ls -lh**
drwxrwxr-x 3 user group 4.0K Sep 18 17:46 cromwell-executions
drwxrwxrwx 2 user group 4.0K Sep 18 19:39 cromwell-workflow-logs
[...snip...]
drwxr-xr-x 4 user group 4.0K Sep 18 15:45 input
        
```

The pipeline outputs can be found in `cromwell-executions` with an
idiosyncratic naming scheme. This directory contains a lot of hard links, so
links have to be preserved when copying back to /data or the size of the folder
will increase substantially.



```


[user@cn3144]$ **cp -ra cromwell-executions $wd**

[user@cn3144]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf]$

```





```


[user@biowulf]$ **sinteractive --cpus-per-task=8 --mem=20g --gres=lscratch:30**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144]$ **wd=$PWD**  # so we can copy results back later
[user@cn3144]$ **cd /lscratch/$SLURM\_JOB\_ID**
[user@cn3144]$ **module load encode-atac-seq-pipeline/1.1.7**
[user@cn3144]$ **cp -Lr ${EASP\_TEST\_DATA:-none}/\* .**
[user@cn3144]$ **tree**
.
|-- [user    738]  ENCSR889WQX.json.1.0
|-- [user    870]  ENCSR889WQX.json.1.1.7
`-- [user   4.0K]  input
    |-- [user   4.0K]  rep1
    |   |-- [user   213M]  ENCFF439VSY_5M.fastq.gz
    |   `-- [user   213M]  ENCFF683IQS_5M.fastq.gz
    `-- [user   4.0K]  rep2
        |-- [user   210M]  ENCFF463QCX_5M.fastq.gz
        `-- [user   211M]  ENCFF992TSA_5M.fastq.gz

[user@cn3144]$ **cat ENCSR889WQX.json.1.1.7**
{
    "atac.pipeline_type" : "atac",
    "atac.genome_tsv" : "/fdb/encode-atac-seq-pipeline/v1/mm10/mm10.tsv",
    "atac.fastqs_rep1_R1" :
        ["input/rep1/ENCFF683IQS_5M.fastq.gz", "input/rep1/ENCFF439VSY_5M.fastq.gz"],
    "atac.fastqs_rep2_R1" :
        ["input/rep2/ENCFF992TSA_5M.fastq.gz", "input/rep2/ENCFF463QCX_5M.fastq.gz"],
    "atac.paired_end" : false,
    "atac.multimapping" : 4,
    "atac.trim_adapter.auto_detect_adapter" : true,
    "atac.smooth_win" : 73,
    "atac.enable_idr" : true,
    "atac.idr_thresh" : 0.05,
    "atac.qc_report.name" : "ENCSR889WQX (subsampled to 5M reads)",
    "atac.qc_report.desc" : "ATAC-seq on Mus musculus C57BL/6 frontal cortex adult",
    "atac.bowtie2_cpu" : 4
}
        
```

In this example the pipeline will only be run locally - i.e. it will not submit
 tasks as slurm jobs



```

[user@cn3144]$ **java -jar -Dconfig.file=$EASP\_BACKEND\_CONF \
 -Dbackend.default=singularity \
 -Dbackend.providers.singularity.config.concurrent-job-limit=2 \
 ${CROMWELL\_JAR} run $EASP\_WDL -i ENCSR889WQX.json.1.1.7 -o $EASP\_WFOPTS \
 -m meta.json** 
[...much output...]
[user@cn3144]$ **ls -lh**
drwxrwxr-x 3 user group 4.0K Sep 18 17:46 cromwell-executions
drwxrwxrwx 2 user group 4.0K Sep 18 19:39 cromwell-workflow-logs
[...snip...]
drwxr-xr-x 4 user group 4.0K Sep 18 15:45 input
-rw-r--r-- 1 user group 283K Sep 18 14:15 meta.json
        
```

The pipeline outputs can be found in `cromwell-executions` with an
idiosyncratic naming scheme. This directory contains a lot of hard links, so
links have to be preserved when copying back to /data or the size of the folder
will increase substantially.



```


[user@cn3144]$ **cp -ra cromwell-executions $wd**

[user@cn3144]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf]$

```





```


[user@biowulf]$ **sinteractive --cpus-per-task=8 --mem=20g --gres=lscratch:30**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144]$ **wd=$PWD**  # so we can copy results back later
[user@cn3144]$ **cd /lscratch/$SLURM\_JOB\_ID**
[user@cn3144]$ **module load encode-atac-seq-pipeline/1.6.1**
[user@cn3144]$ **cp -Lr ${EASP\_TEST\_DATA:-none}/\* .**
[user@cn3144]$ **tree**
.
|-- [user    738]  ENCSR889WQX.json.1.0
|-- [user    870]  ENCSR889WQX.json.1.1.7
|-- [user    870]  ENCSR889WQX.json.1.6.1
|-- [user    870]  ENCSR889WQX.json.1.9.0
`-- [user   4.0K]  input
    |-- [user   4.0K]  rep1
    |   |-- [user   213M]  ENCFF439VSY_5M.fastq.gz
    |   `-- [user   213M]  ENCFF683IQS_5M.fastq.gz
    `-- [user   4.0K]  rep2
        |-- [user   210M]  ENCFF463QCX_5M.fastq.gz
        `-- [user   211M]  ENCFF992TSA_5M.fastq.gz

[user@cn3144]$ **cat ENCSR889WQX.json.1.6.1**
{
    "atac.pipeline_type" : "atac",
    "atac.genome_tsv" : "/fdb/encode-atac-seq-pipeline/v1/mm10/mm10.tsv",
    "atac.mito_chr_name": "chrM",
    "atac.fastqs_rep1_R1" :
        ["input/rep1/ENCFF683IQS_5M.fastq.gz", "input/rep1/ENCFF439VSY_5M.fastq.gz"],
    "atac.fastqs_rep2_R1" :
        ["input/rep2/ENCFF992TSA_5M.fastq.gz", "input/rep2/ENCFF463QCX_5M.fastq.gz"],
    "atac.paired_end" : false,
    "atac.multimapping" : 4,
    "atac.auto_detect_adapter" : true,
    "atac.smooth_win" : 73,
    "atac.enable_idr" : true,
    "atac.idr_thresh" : 0.05,
    "atac.title" : "ENCSR889WQX (subsampled to 5M reads)",
    "atac.description" : "ATAC-seq on Mus musculus C57BL/6 frontal cortex adult"
}
        
```

In this example the pipeline will only be run locally - i.e. it will not submit
 tasks as slurm jobs. Follow the [caper](https://github.com/ENCODE-DCC/caper)
 docs to set up a config file for slurm submission. This has to be done only once.



```

[user@cn3144]$ **mkdir -p ~/.caper && caper init local**
[user@cn3144]$ **caper run $EASP\_WDL -i ENCSR889WQX.json.1.6.1**
[...much output...]
This workflow ran successfully. There is nothing to troubleshoot

```

This version of the pipeline comes with a tool to copy and organize pipeline output



```

[user@cn3144]$ **croo --method copy --out-dir=${wd}/ENCSR889WQX \
 atac/18a97503-94e1-4d75-a1e6-6a582a4c5407/metadata.json**

```




Note the switch to annotation v3.



```

[user@biowulf]$ **sinteractive --cpus-per-task=8 --mem=20g --gres=lscratch:30**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144]$ **wd=$PWD**  # so we can copy results back later
[user@cn3144]$ **cd /lscratch/$SLURM\_JOB\_ID**
[user@cn3144]$ **module load encode-atac-seq-pipeline/1.9.0**
[user@cn3144]$ **cp -Lr ${EASP\_TEST\_DATA:-none}/\* .**
[user@cn3144]$ **tree**
.
|-- [user    738]  ENCSR889WQX.json.1.0
|-- [user    870]  ENCSR889WQX.json.1.1.7
|-- [user    870]  ENCSR889WQX.json.1.6.1
|-- [user    870]  ENCSR889WQX.json.1.9.0
`-- [user   4.0K]  input
    |-- [user   4.0K]  rep1
    |   |-- [user   213M]  ENCFF439VSY_5M.fastq.gz
    |   `-- [user   213M]  ENCFF683IQS_5M.fastq.gz
    `-- [user   4.0K]  rep2
        |-- [user   210M]  ENCFF463QCX_5M.fastq.gz
        `-- [user   211M]  ENCFF992TSA_5M.fastq.gz

[user@cn3144]$ **cat ENCSR889WQX.json.1.9.0**
{
    "atac.pipeline_type" : "atac",
    "atac.genome_tsv" : "/fdb/encode-atac-seq-pipeline/v3/mm10/mm10.tsv",
    "atac.mito_chr_name": "chrM",
    "atac.fastqs_rep1_R1" :
        ["input/rep1/ENCFF683IQS_5M.fastq.gz", "input/rep1/ENCFF439VSY_5M.fastq.gz"],
    "atac.fastqs_rep2_R1" :
        ["input/rep2/ENCFF992TSA_5M.fastq.gz", "input/rep2/ENCFF463QCX_5M.fastq.gz"],
    "atac.paired_end" : false,
    "atac.multimapping" : 4,
    "atac.auto_detect_adapter" : true,
    "atac.smooth_win" : 73,
    "atac.enable_idr" : true,
    "atac.idr_thresh" : 0.05,
    "atac.title" : "ENCSR889WQX (subsampled to 5M reads)",
    "atac.description" : "ATAC-seq on Mus musculus C57BL/6 frontal cortex adult"
}
        
```

In this example the pipeline will only be run locally - i.e. it will not submit
 tasks as slurm jobs. If you wish to run in Slurm execution mode, you will need to create
 a custom Slurm configuration following the [caper](https://github.com/ENCODE-DCC/caper)
 docs to set up a config file for slurm submission. We recommend sticking to local execution
 unless Slurm submission becomes necessary.



```

[user@cn3144]$ **mkdir -p ~/.caper && caper init local**
[user@cn3144]$ **caper run $EASP\_WDL -i ENCSR889WQX.json.1.9.0**
[...much output...]
This workflow ran successfully. There is nothing to troubleshoot

```

This version of the pipeline comes with a tool to copy and organize pipeline output.



```

[user@cn3144]$ **ls atac**
a0fb9f58-ede3-4c02-9bcc-26d21ab5ccbb
[user@cn3144]$ **croo --method copy --out-dir=${wd}/ENCSR889WQX \
 atac/a0fb9f58-ede3-4c02-9bcc-26d21ab5ccbb/metadata.json**

```




Continues to use the v3 annotation. However, caper apparently changed
 significantly so you should backup your old caper configuration and create
 fresh config files for this version.



```

[user@biowulf]$ **sinteractive --cpus-per-task=8 --mem=20g --gres=lscratch:30**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144]$ **wd=$PWD**  # so we can copy results back later
[user@cn3144]$ **cd /lscratch/$SLURM\_JOB\_ID**
[user@cn3144]$ **module load encode-atac-seq-pipeline/2.1.0**
[user@cn3144]$ **cp -Lr ${EASP\_TEST\_DATA:-none}/\* .**
[user@cn3144]$ **tree**
.
|-- [user    738]  ENCSR889WQX.json.1.0
|-- [user    870]  ENCSR889WQX.json.1.1.7
|-- [user    870]  ENCSR889WQX.json.1.6.1
|-- [user    870]  ENCSR889WQX.json.1.9.0
|-- [user    870]  ENCSR889WQX.json.2.1.0
`-- [user   4.0K]  input
    |-- [user   4.0K]  rep1
    |   |-- [user   213M]  ENCFF439VSY_5M.fastq.gz
    |   `-- [user   213M]  ENCFF683IQS_5M.fastq.gz
    `-- [user   4.0K]  rep2
        |-- [user   210M]  ENCFF463QCX_5M.fastq.gz
        `-- [user   211M]  ENCFF992TSA_5M.fastq.gz

[user@cn3144]$ **cat ENCSR889WQX.json.2.1.0**
{
    "atac.pipeline_type" : "atac",
    "atac.genome_tsv" : "/fdb/encode-atac-seq-pipeline/v3/mm10/mm10.tsv",
    "atac.mito_chr_name": "chrM",
    "atac.fastqs_rep1_R1" :
        ["input/rep1/ENCFF683IQS_5M.fastq.gz", "input/rep1/ENCFF439VSY_5M.fastq.gz"],
    "atac.fastqs_rep2_R1" :
        ["input/rep2/ENCFF992TSA_5M.fastq.gz", "input/rep2/ENCFF463QCX_5M.fastq.gz"],
    "atac.paired_end" : false,
    "atac.multimapping" : 4,
    "atac.auto_detect_adapter" : true,
    "atac.smooth_win" : 73,
    "atac.enable_idr" : true,
    "atac.idr_thresh" : 0.05,
    "atac.title" : "ENCSR889WQX (subsampled to 5M reads)",
    "atac.description" : "ATAC-seq on Mus musculus C57BL/6 frontal cortex adult"
}
        
```

In this example the pipeline will only be run locally - i.e. it will not submit
 tasks as slurm jobs. Follow the [caper](https://github.com/ENCODE-DCC/caper)
 docs to set up a config file for slurm submission. This has to be done only once.



```

[user@cn3144]$ **[[ -d ~/.caper ]] && mv ~/.caper ~/caper.$(date +%F).bak** # back up old caper config
[user@cn3144]$ **mkdir -p ~/.caper && caper init local**
[user@cn3144]$ # note the need for --singularity in this version
[user@cn3144]$ **caper run $EASP\_WDL -i ENCSR889WQX.json.2.1.0 --singularity**
[...much output...]
This workflow ran successfully. There is nothing to troubleshoot

```

This version of the pipeline comes with a tool to copy and organize pipeline output.



```

[user@cn3144]$ **ls atac**
a0fb9f58-ede3-4c02-9bcc-26d21ab5ccbb
[user@cn3144]$ **croo --method copy --out-dir=${wd}/ENCSR889WQX \
 atac/a0fb9f58-ede3-4c02-9bcc-26d21ab5ccbb/metadata.json**

```


 

Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file. For example the following batch job will run a local job (assuming the caper config file is set up correctly):



```

#! /bin/bash

wd=$PWD
module load encode-atac-seq-pipeline/2.2.0 || exit 1

cd /lscratch/$SLURM_JOB_ID

mkdir input
cp -rL $EASP_TEST_DATA/* .
caper run $EASP_WDL -i ENCSR889WQX.json.2.1.0
rc=$?
croo --method copy --out-dir=${wd}/ENCSR889WQX \
    atac/*/metadata.json
exit $rc

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --time=4:00:00 --cpus-per-task=8 --mem=20g --gres=lscratch:50 encode-atac-seq-pipeline.sh
```







