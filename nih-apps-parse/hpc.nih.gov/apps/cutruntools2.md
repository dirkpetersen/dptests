














cutruntools2 on Biowulf


  







|  |  |  |
| --- | --- | --- |
|  | 
[Biowulf High Performance Computing at the NIH](https://hpc.nih.gov)
 | 





[GitHub](https://github.com/NIH-HPC)
[YouTube](https://www.youtube.com/channel/UCx-kNd1kBskYr5KLT9-Erew)
[@nih_hpc](https://twitter.com/nih_hpc)
[RSS Feed](https://hpc.nih.gov/hpc_RSS.xml)
 |




* [Status](https://hpc.nih.gov/systems/status/)
* [Applications](https://hpc.nih.gov/apps/)
* [Reference Data](https://hpc.nih.gov/apps/db.php)
* [Storage](https://hpc.nih.gov/storage/)
* [User Guides](https://hpc.nih.gov/docs/user_guides.html)
* [Training](https://hpc.nih.gov/training/)
* [User Dashboard](https://hpc.nih.gov/dashboard/)
* [How To](https://hpc.nih.gov/docs/how_to.html)
* [About](https://hpc.nih.gov/about/)







cutruntools2 on Biowulf


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



cutruntools2 is a major update of [CutRunTools](https://hpc.nih.gov/apps/CutRunTools.html), including a set of new features specially designed for CUT&RUN and CUT&Tag experiments. Both of the bulk and single-cell data can be processed, analyzed and interpreted. 






### References:


* Fulong Yu, Vijay G Sankaran, Guo-Cheng Yuan*CUT&RUNTools 2.0:a pipeline for single-cell and bulk-level CUT&RUN and CUT&Tag data analysis* Bioinformatics, 2021
 [PubMed](https://pubmed.ncbi.nlm.nih.gov/34244724/) | 
 [Journal](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btab507/6318389)


Documentation
* cutruntools2 Github:[Github](https://github.com/fl-yu/CUT-RUNTools-2.0)


Important Notes
* Module Name: cutruntools2 (see [the modules page](/apps/modules.html) for more information)
 * cutruntools2 needs config JSON file as input, please do not edit the software\_config, but only modify input\_output, and motif\_finding to match your username with working path according. This is an example for bulk-config.json:
 
```

{
    "software_config": {
        "Rscriptbin": "/usr/bin",
        "pythonbin": "/opt/conda/envs/app/bin",
        "perlbin": "/opt/conda/envs/app/bin",
        "javabin": "/usr/bin",
        "bowtie2bin": "/opt/conda/envs/app/bin",
        "samtoolsbin": "/opt/conda/envs/app/bin",
        "macs2bin": "/opt/conda/envs/app/bin",
        "memebin": "/opt/conda/envs/meme/bin",
        "bedopsbin": "/opt/conda/envs/app/bin",
        "bedtoolsbin": "/opt/conda/envs/app/bin",
        "path_deeptools": "/opt/conda/envs/app/bin",
        "path_parallel": "/opt/conda/envs/app/bin",
        "path_tabix": "/opt/conda/envs/app/bin",
        "bt2idx": "/fdb/igenomes/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index",
        "genome_sequence": "/fdb/igenomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa",
        "extratoolsbin": "/opt/CUT-RUNTools-2.0/install",
        "extrasettings": "/opt/CUT-RUNTools-2.0/install",
        "kseqbin": "/opt/CUT-RUNTools-2.0/install",
        "adapterpath": "/opt/CUT-RUNTools-2.0/adapters",
        "trimmomaticbin": "/opt/CUT-RUNTools-2.0/install",
        "picardbin": "/opt/CUT-RUNTools-2.0/install",
        "picardjarfile": "picard-2.8.0.jar",
        "trimmomaticjarfile": "trimmomatic-0.36.jar",
	"makecutmatrixbin": "/opt/conda/envs/app/bin"
    },
    "input_output": {
        "fastq_directory": "/data/apptest1/cutruntools2/exampleData",
        "workdir": "/data/apptest1/cutruntools2/bulk-example-test",
        "fastq_sequence_length": "42",
        "organism_build": "hg38",
        "spike_in": "FALSE",
        "spike_in_norm": "FALSE",
        "spikein_scale": "10000",
        "frag_120": "TRUE",
        "peak_caller": "macs2",
        "dup_peak_calling": "FALSE",
        "cores": "8",
        "experiment_type": "CUT&RUN"
    },
    "motif_finding": {
        "num_bp_from_summit": "100",
        "num_peaks": "1000",
        "total_peaks": "2000",
        "motif_scanning_pval": "0.0005",
        "num_motifs": "10"
    }
}


	
```




Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive --ntasks=20 --cpus-per-task=1 --mem=20G --nodes=1**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load cutruntools2**
[user@cn3144 ~]$ **mkdir /data/$USER/cutruntools2/**
[user@cn3144 ~]$ **cd /data/$USER/cutruntools2/**
[user@cn3144 ~]$ **cp -r ${CUTRUNTOOLS\_TEST\_DATA:-none}/\* .**
[user@cn3144 ~]$ **sed -i "s/apptest1/$USER/g" \*config.json** # replace username with yours
[user@cn3144 ~]$ **python ./validate.py bulk-config.json --ignore-input-output --software**
[user@cn3144 ~]$ **bash run\_bulkModule.sh bulk-config.json test1 &>test1\_bulk.out**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. cutruntools2\_sc.sh). For example:



```


#!/bin/bash
module load cutruntools2
mkdir -p /data/$USER/cutruntools2
cd /data/$USER/cutruntools2
cp -r ${CUTRUNTOOLS_TEST_DATA:-none}/* .
sed -i "s/apptest1/$USER/g" *config.json # replace username with yours
python ./validate.py sc-config.json --ignore-input-output --software
bash run_scModule.sh sc-config.json


```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --ntasks=20 --cpus-per-task=1 --mem=70G --nodes=1 cutruntools2_sc.sh
```












[HPC @ NIH](https://hpc.nih.gov/)  ~
[Contact](/about/contact.html)


[Disclaimer](/docs/disclaimer.html) ~
[Privacy](/docs/privacy.html) ~
[Accessibility](/docs/accessibility.html) ~
[CIT](http://cit.nih.gov/) ~
[NIH](http://www.nih.gov/) ~
[DHHS](http://www.dhhs.gov/) ~
[USA.gov](http://www.firstgov.gov/)



  


