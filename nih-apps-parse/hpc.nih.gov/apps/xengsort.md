

document.querySelector('title').textContent = 'xengsort on Biowulf';
xengsort on Biowulf


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



xengsort is a fast xenograft read sorter based on space-efficient k-mer hashing.



### References:


* [Zentgraf, Jens, and Sven Rahmann. "Fast lightweight accurate xenograft sorting." *bioRxiv* (2020).](https://www.biorxiv.org/content/biorxiv/early/2020/05/19/2020.05.14.095604.full.pdf)


Documentation
* [xengsort on GitLab](https://gitlab.com/genomeinformatics/xengsort/)


Important Notes
* Module Name: xengsort (see [the modules page](/apps/modules.html) for more information)
 * This application is designed to work with [snakemake](https://hpc.nih.gov/apps/snakemake.html).



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf ~]$ **sinteractive --cpus-per-task=8 --mem=32g**
salloc.exe: Pending job allocation 61524097
salloc.exe: job 61524097 queued and waiting for resources
salloc.exe: job 61524097 has been allocated resources
salloc.exe: Granted job allocation 61524097
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3137 are ready for job
srun: error: x11: no local DISPLAY defined, skipping

[user@cn3137 ~]$ **cd /data/$USER**

[user@cn3137 user]$ **git clone https://gitlab.com/genomeinformatics/xengsort.git**
Cloning into 'xengsort'...
remote: Enumerating objects: 226, done.
remote: Counting objects: 100% (226/226), done.
remote: Compressing objects: 100% (103/103), done.
remote: Total 226 (delta 139), reused 183 (delta 117), pack-reused 0
Receiving objects: 100% (226/226), 108.80 KiB | 0 bytes/s, done.
Resolving deltas: 100% (139/139), done.

[user@cn3137 user]$ **cd xengsort/**

[user@cn3137 xengsort]$ **module load xengsort snakemake**
[+] Loading xengsort  28762aac  on cn3137
[+] Loading singularity  3.5.3  on cn3137
[+] Loading snakemake  5.19.3

[user@cn3137 xengsort]$ **snakemake -j 8**
Building DAG of jobs...
Using shell: /usr/bin/bash
Provided cores: 8
Rules claiming more threads will be scaled down.
Job counts:
        count   jobs
        1       all
        1       build_index
        1       classify_mouse_exomes
        2       download_mouse_exomes
        1       download_refs
        6

[Mon Jul 20 14:15:43 2020]
rule download_refs:
    output: ref/Homo_sapiens.GRCh38.dna.toplevel.fa.gz, ref/Mus_musculus.GRCm38.dna.toplevel.fa.gz, ref/Homo_sapiens.GRCh38.cdna.all.fa.gz, ref/Mus_musculus.GRCm38.cdna.all.fa.gz
    jobid: 5


[Mon Jul 20 14:15:43 2020]
rule download_mouse_exomes:
    output: raw/BALBc-M1-normal_1.fq.gz.1
    jobid: 3
    wildcards: filename=BALBc-M1-normal_1.fq.gz.1


[Mon Jul 20 14:15:43 2020]
rule download_mouse_exomes:
    output: raw/BALBc-M1-normal_2.fq.gz.1
    jobid: 4
    wildcards: filename=BALBc-M1-normal_2.fq.gz.1

--2020-07-20 14:15:43--  https://sra-pub-src-1.s3.amazonaws.com/SRR9130497/BALBc-M1-normal_2.fq.gz.1
--2020-07-20 14:15:43--  https://sra-pub-src-1.s3.amazonaws.com/SRR9130497/BALBc-M1-normal_1.fq.gz.1
--2020-07-20 14:15:43--  ftp://ftp.ensembl.org/pub/release-98/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
Resolving dtn05-e0 (dtn05-e0)... 10.1.200.241
Connecting to dtn05-e0 (dtn05-e0)|10.1.200.241|:3128... connected.
Proxy request sent, awaiting response... Resolving dtn05-e0 (dtn05-e0)... Resolving dtn05-e0 (dtn05-e0)... 10.1.200.24110.1.200.241
Connecting to dtn05-e0 (dtn05-e0)|10.1.200.241|:3128...
Connecting to dtn05-e0 (dtn05-e0)|10.1.200.241|:3128... connected.
connected.
Proxy request sent, awaiting response... Proxy request sent, awaiting response... 200 OK
Length: 5034316763 (4.7G) [application/x-troff-man]
Saving to: ‘raw/BALBc-M1-normal_2.fq.gz.1’

 0% [                                                             ] 0           --.-K/s              200 OK
Length: 4783409755 (4.5G) [application/x-troff-man]
Saving to: ‘raw/BALBc-M1-normal_1.fq.gz.1’

 0% [                                                             ] 33,550,787  32.2MB/s             200 Gatewaying
Length: 1107654500 (1.0G) [text/plain]
Saving to: ‘ref/Homo_sapiens.GRCh38.dna.toplevel.fa.gz’
[...snip...]

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. xengsort.sh). For example:



```

#!/bin/bash
set -e
cd /data/${USER}
git clone https://gitlab.com/genomeinformatics/xengsort.git
cd xengsort/
module load xengsort snakemake
snakemake -j 8

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] xengsort.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. xengsort.swarm). For example:



```

cd /path/to/snakefile1 && snakemake -j 8
cd /path/to/snakefile2 && snakemake -j 8
cd /path/to/snakefile3 && snakemake -j 8
cd /path/to/snakefile4 && snakemake -j 8

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f xengsort.swarm [-g #] [-t #] --module xengsort snakemake
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module xengsort Loads the xengsort module for each subjob in the swarm 
 | |
 | |
 | |








