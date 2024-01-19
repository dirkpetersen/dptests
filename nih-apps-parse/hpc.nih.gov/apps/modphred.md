

document.querySelector('title').textContent = 'Modphred on Biowulf';
Modphred on Biowulf


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


ModPhred is a pipeline for detection, annotation and visualisation of DNA/RNA modifications (From Authors' documentation).. 


### References:


* Pryszcz LP and Novoa EM
*[ModPhred: an integrative toolkit for the analysis and storage of nanopore sequencing DNA and RNA modification data](https://pubmed.ncbi.nlm.nih.gov/34293115/)*. Bioinformatics, 38:257-260 (2022)


Documentation
* Modphred Main Site: [ReadTheDocs](https://modphred.readthedocs.io/en/latest/index.html)


Important Notes
* Module Name: modphred (see [the modules page](/apps/modules.html) for more information)
* You will need to request GPU resources to run modphred (see example below). the current version of modphred will not work on a100 GPUs



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --gres=gpu:p100:1,lscratch:200 --mem=16g --cpus-per-task=6**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn4224 are ready for job

[user@cn4224 ~]$ **export TMPDIR=/lscratch/${SLURM\_JOB\_ID}**
[user@cn4224 ~]$ **module load modphred**
[user@cn4224 ~]$ **cd /data/${USER}**
[user@cn4224 ~]$ **wget https://public-docs.crg.es/enovoa/public/lpryszcz/src/modPhred/test/ -q -r -c -nc -np -nH --cut-dirs=6 --reject="index.html\*"**
[user@cn4224 ~]$ **run -f ref/ECOLI.fa -o OUTPUT -i PRJEB22772/\* -t4 --host /usr/bin/guppy\_basecall\_server**
[2022-12-29 11:34:36] ===== Welcome, welcome to modPhred pipeline! =====
[2022-12-29 11:34:36] Starting /usr/bin/guppy_basecall_server ... [mem:   134 MB]
[2022-12-29 11:34:40] Encoding modification info from 2 directories... [mem:   134 MB]
[2022-12-29 11:34:40]  PRJEB22772/MARC_ZFscreens_R9.4_1D-Ecoli-run_FAF05145 with 4 Fast5 file(s)... [mem:   134 MB]
[2022-12-29 11:36:21]   DNA alphabet with 2 modification(s) {'A': ['Y'], 'C': ['Z'], 'G': [], 'T': []}. symbol2modbase: {'Y': '6mA', 'Z': '5mC'}
[2022-12-29 11:39:20]   106,722,317 bases saved in FastQ, of those: 332,925 6mA [ 0.312%], 160,288 5mC [ 0.150%]   
[2022-12-29 11:39:20]  PRJEB22772/MARC_ZFscreens_R9.4_2D-Ecoli-run_FAF05711 with 1 Fast5 file(s)... [mem:   988 MB]
[2022-12-29 11:41:40]   DNA alphabet with 2 modification(s) {'A': ['Y'], 'C': ['Z'], 'G': [], 'T': []}. symbol2modbase: {'Y': '6mA', 'Z': '5mC'}
[2022-12-29 11:41:40]   29,420,071 bases saved in FastQ, of those: 91,315 6mA [ 0.310%], 51,017 5mC [ 0.173%]   
[2022-12-29 11:41:40] Aligning FastQ files from 2 directories... [mem:  2292 MB]
[2022-12-29 11:41:40]   > OUTPUT/minimap2/MARC_ZFscreens_R9.4_1D-Ecoli-run_FAF05145.bam [mem:  2292 MB]
[2022-12-29 11:42:01]   > OUTPUT/minimap2/MARC_ZFscreens_R9.4_2D-Ecoli-run_FAF05711.bam [mem:  2292 MB]
[2022-12-29 11:42:09] Indexing bam file(s)... [mem:  2392 MB]
[2022-12-29 11:42:11] Reporting positions that are likely modified to OUTPUT/mod.gz ... [mem:  2392 MB]
[2022-12-29 11:42:11]  Getting regions covered by at least 25 reads... [mem:  2392 MB]
[2022-12-29 11:42:12]   7 regions to process... [mem:  2392 MB]
[2022-12-29 11:45:50] Loading modification data... [mem:  2392 MB]
[2022-12-29 11:45:50] Plotting... [mem:  2392 MB]
[2022-12-29 11:45:51] Saving modified positions with max frequency as OUTPUT/mod.bed (bedMethyl file) ... [mem:  2392 MB]
[2022-12-29 11:45:51]  and separately for every BAM file as OUTPUT/minimap2/*.bed ... [mem:  2392 MB]
[2022-12-29 11:46:00] Saving plots for depth, basecall_accuracy, mod_frequency, median_mod_prob to OUTPUT/plots ... [mem:  2392 MB]
[2022-12-29 11:46:05] You can remove reads directory: rm -r OUTPUT/reads/ [mem:  2392 MB]
[2022-12-29 11:46:05] All finished! Have a nice day :) [mem:  2392 MB]
#Time elapsed: 0:11:28.435555

[user@cn4224 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. modphred.sh) similar to the following.



```

#! /bin/bash

module load modphred
run -f ref/ECOLI.fa -o OUTPUT -i PRJEB22772/* -t4 --host /usr/bin/guppy_basecall_server


```

Submit these jobs using the Slurm [sbatch](/docs/userguide.html) command:



```
sbatch --partition=gpu --cpus-per-task=6 --mem=16g --gres=lscratch:200,gpu:p100:1 modphred.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile for the pipeline (e.g. modphred.swarm). For example:



```

run -f ref/ECOLI.fa -o OUTPUT2 -i PRJEB22772/* -t4 --host /usr/bin/guppy_basecall_server
run -f ref/ECOLI.fa -o OUTPUT3 -i PRJEB22773/* -t4 --host /usr/bin/guppy_basecall_server
run -f ref/ECOLI.fa -o OUTPUT4 -i PRJEB22774/* -t4 --host /usr/bin/guppy_basecall_server

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f modphred.swarm --partition=gpu -g 16 -t 6 --gres=gpu:p100:1,lscratch:200 --module modphred
```

where


|  |  |  |  |
| --- | --- | --- | --- |
| -g #  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | --module modphred  Loads the modphred module for each subjob in the swarm
 | |
 | |








