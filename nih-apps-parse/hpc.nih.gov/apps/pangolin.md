

document.querySelector('title').textContent = 'pangolin on Biowulf';
pangolin on Biowulf


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



Phylogenetic Assignment of Named Global Outbreak LINeages. pangolin was developed to implement the dynamic nomenclature of SARS-CoV-2 lineages, known as the Pango nomenclature. It allows a user to assign a SARS-CoV-2 genome sequence the most likely lineage (Pango lineage) to SARS-CoV-2 query sequences.



Documentation
* [pangolin Main Site](https://cov-lineages.org/pangolin.html)
* [pangolin on GitHub](https://github.com/cov-lineages/pangolin)


Important Notes
* Module Name: pangolin (see [the modules page](/apps/modules.html) for more information)
 * pangolin should be run with lscratch and the --tempdir option should be used to make sure that temporary data is written to the correct location. See example below.
 * Environment variables set 
	+ PATH
	+ PANGOLIN\_TESTDATA



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf ~]$ **sinteractive -c4 --mem=8g --gres=lscratch:10**
salloc.exe: Pending job allocation 12440727
salloc.exe: job 12440727 queued and waiting for resources
salloc.exe: job 12440727 has been allocated resources
salloc.exe: Granted job allocation 12440727
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn0913 are ready for job
srun: error: x11: no local DISPLAY defined, skipping
error: unable to open file /tmp/slurm-spank-x11.12440727.0
slurmstepd: error: x11: unable to read DISPLAY value

[user@cn0913 ~]$ **cd /lscratch/$SLURM\_JOB\_ID**

[user@cn0913 12440727]$ **module load pangolin**
[+] Loading pangolin  2.3.6  on cn0913
[+] Loading singularity  3.7.3  on cn0913

[user@cn0913 12440727]$ **cp $PANGOLIN\_TESTDATA/\* .**

[user@cn0913 12440727]$ **mkdir tempdir**

[user@cn0913 12440727]$ **pangolin --tempdir=./tempdir cluster.fasta**
Found the snakefile
The query file is:/lscratch/12440727/cluster.fasta
EDB003  sequence too short
EDB004  has an N content of 0.98
Looking in /opt/conda/envs/pangolin/lib/python3.7/site-packages/pangoLEARN/data for data files...

Data files found
Trained model:  /opt/conda/envs/pangolin/lib/python3.7/site-packages/pangoLEARN/data/decisionTree_v1.joblib
Header file:    /opt/conda/envs/pangolin/lib/python3.7/site-packages/pangoLEARN/data/decisionTreeHeaders_v1.joblib
Lineages csv:   /opt/conda/envs/pangolin/lib/python3.7/site-packages/pangoLEARN/data/lineages.metadata.csv
Job counts:
        count   jobs
        1       add_failed_seqs
        1       align_to_reference
        1       all
        1       minimap2_check_distance
        1       overwrite
        1       pangolearn
        1       parse_paf
        1       type_variants_b117
        1       type_variants_b12142
        1       type_variants_b1351
        1       type_variants_p1
        1       type_variants_p2
        1       type_variants_p3
        13
Job counts:
        count   jobs
        1       parse_paf
        1
[M::mm_idx_gen::0.003*1.49] collected minimizers
[M::mm_idx_gen::0.004*1.33] sorted minimizers
[M::main::0.004*1.33] loaded/built the index for 1 target sequence(s)
[M::mm_mapopt_update::0.004*1.30] mid_occ = 100
[M::mm_idx_stat] kmer size: 19; skip: 19; is_hpc: 0; #seq: 1
[M::mm_idx_stat::0.004*1.28] distinct minimizers: 2952 (100.00% are singletons); average occurrences: 1.000; average spacing: 10.130
warning: using --pad without --trim has no effect
[M::worker_pipeline::0.095*1.00] mapped 8 sequences
[M::main] Version: 2.17-r941
[M::main] CMD: minimap2 -a -x asm5 -t 1 /opt/conda/envs/pangolin/lib/python3.7/site-packages/pangolin-2.3.6-py3.7.egg/pangolin/data/reference.fasta /lscratch/12440727/./tempdir/tmpat6kcm3p/mappable.fasta
[M::main] Real time: 0.096 sec; CPU: 0.096 sec; Peak RSS: 0.003 GB
loading model 04/09/2021, 13:23:32
processing block of 8 sequences 04/09/2021, 13:23:33
complete 04/09/2021, 13:23:34
Job counts:
        count   jobs
        1       add_failed_seqs
        1
Job counts:
        count   jobs
        1       overwrite
        1
Output file written to: /lscratch/12440727/lineage_report.csv

[user@cn0913 12440727]$ **cat lineage\_report.csv**
taxon,lineage,probability,pangoLEARN_version,status,note
EDB001,B.1.1.26,1.0,2021-03-29,passed_qc,
EDB002,B.1.1.26,1.0,2021-03-29,passed_qc,
EDB005,B,1.0,2021-03-29,passed_qc,
EDB006,B.1.1.26,1.0,2021-03-29,passed_qc,
EDB007,B.1.1.26,1.0,2021-03-29,passed_qc,
EDB008,B.1.1.26,1.0,2021-03-29,passed_qc,
EDB009,B.1.1.26,1.0,2021-03-29,passed_qc,
EDB010,B,1.0,2021-03-29,passed_qc,
EDB003,None,0,2021-03-29,fail,seq_len:2997
EDB004,None,0,2021-03-29,fail,N_content:0.98

[user@cn0913 12440727]$ **exit**
exit
salloc.exe: Relinquishing job allocation 12440727

[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. pangolin.sh). For example:



```

#!/bin/bash
set -e
module load pangolin
pagolin --tempdir=/lscratch/$SLURM_JOB_ID my.fasta

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] [--gres=lscratch:#] pangolin.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. pangolin.swarm). For example:



```

pagolin --tempdir=/lscratch/$SLURM_JOB_ID A.fasta
pagolin --tempdir=/lscratch/$SLURM_JOB_ID B.fasta
pagolin --tempdir=/lscratch/$SLURM_JOB_ID C.fasta
pagolin --tempdir=/lscratch/$SLURM_JOB_ID D.fasta

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f pangolin.swarm [-g #] [-t #] [--gres=lscratch:#] --module pangolin
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module pangolin Loads the pangolin module for each subjob in the swarm 
 | |
 | |
 | |








