

document.querySelector('title').textContent = 'distiller-nf on Biowulf';
distiller-nf on Biowulf


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



A modular Hi-C mapping pipeline for reproducible data analysis, it was used for Micro-C analysis too.
The distiller pipeline aims to provide the following functionality:

- Align the sequences of Hi-C molecules to the reference genome
- Parse .sam alignment and form files with Hi-C pairs
- Filter PCR duplicates
- Aggregate pairs into binned matrices of Hi-C interactions






### References:


* Krietenstein N, Abraham S, Venev SV, Abdennur N, Gibcus J, Hsieh TS, Parsi KM, Yang L, Maehr R, Mirny LA, Dekker J, Rando OJ. *Ultrastructural Details of Mammalian Chromosome Architecture.* Mol Cell. 2020 May 7
 [PubMed](https://pubmed.ncbi.nlm.nih.gov/32213324/) | 
 [Journal](https://www.sciencedirect.com/science/article/pii/S1097276520301519?via%3Dihub)


Documentation
* distiller-nf Github:[Github](https://github.com/open2c/distiller-nf)


Important Notes
* Module Name: distiller-nf (see [the modules page](/apps/modules.html) for more information)
 * distiller-nf is a nextflow pipeline, it requires three configuration files: project.yml(for input and ref seq), nextflow.config(already setted up on Biowulf), and local/cluster.config(Biowulf provides templates which you can download and modify).
 
```

	cp -r ${DISTILLER_CONFIG:-none} .
	
```
* Then to give a customized config file, you have to provide the full path to custom.config.
 
```

	-profile custom --custom_config /your_full_path_to/custom.config
	
```
* To run test data with test\_project.yml, you also need to save a local copy of the test data. 
 
```

	cp -r ${DISTILLER_TEST_DATA:-none} .
	
```




Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive --cpus-per-task=10 --mem=10G**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144]$ **module load distiller-nf**
[user@cn3144]$ **mkdir /data/$USER/distiller/**
[user@cn3144]$ **cd /data/$USER/distiller/**
[user@cn3144]$ **cp -r ${DISTILLER\_TEST\_DATA:-none} .**
[user@cn3144]$ **cp -r ${DISTILLER\_CONFIG:-none} .**
[user@cn3144]$ **nextflow $DISTILLER/distiller.nf -params-file ./test/test\_project.yml**
N E X T F L O W  ~  version 21.04.1
Launching `/usr/local/apps/distiller-nf/0.3.3/distiller.nf` [trusting_euler] - revision: a0539e1286
executor >  local (26)
[-        ] process > download_truncate_chunk_fastqs                                      -
[-        ] process > local_truncate_chunk_fastqs                                         -
[-        ] process > fastqc                                                              -
[a8/eb3aac] process > map_parse_sort_chunks (library:MATa_R1 run:lane1 chunk:0)           [100%] 5 of 5 ✔
[0e/e85c28] process > merge_dedup_splitbam (library:MATa_R2)                              [100%] 4 of 4 ✔
[50/2f4035] process > bin_zoom_library_pairs (library:MATa_R2 filter:mapq_30)             [100%] 8 of 8 ✔
[d1/bd0918] process > merge_zoom_library_group_coolers (library_group:all filter:mapq_30) [100%] 6 of 6 ✔
[3a/d421e5] process > merge_stats_libraries_into_groups (library_group:MATalpha)          [100%] 3 of 3 ✔
Completed at: 04-Aug-2021 12:44:36
Duration    : 1m 14s
CPU hours   : 0.5
Succeeded   : 26



[user@cn3144]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. distiller-nf.sh). For example:



```


#!/bin/bash
set -e
module load distiller-nf
cd /data/$USER/distiller/
cp -r ${DISTILLER_TEST_DATA:-none} .
cp -r ${DISTILLER_CONFIG:-none} .
nextflow $DISTILLER/distiller.nf -params-file ./test/test_project.yml -profile cluster


```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch distiller-nf.sh
```









