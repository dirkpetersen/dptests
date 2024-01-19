

document.querySelector('title').textContent = 'verkko on Biowulf';
verkko on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
 |



From the tool description




> 
> Verkko is a hybrid genome assembly pipeline developed for telomere-to-telomere
> assembly of PacBio HiFi and Oxford Nanopore reads. Verkko is Finnish for net,
> mesh and graph.
> 
> Verkko uses Canu to correct remaining errors in the HiFi reads, builds a
> multiplex de Bruijn graph using MBG, aligns the Oxford Nanopore reads to the
> graph using GraphAligner, progressively resolves loops and tangles first with
> the HiFi reads then with the aligned Oxford Nanopore reads, and finally creates
> contig consensus sequences using Canu's consensus module.
> 


### References:


* M. Rautiainen, S. Nurk, B. P. Walenz, G. A. Logsdon, D. Porubsky, A. Rhie, E. E. Eichler, A. M. Phillippy, S. Koren. 
 *Verkko: telomere-to-telomere assembly of diploid chromosomes.*
[bioRxiv](https://doi.org/10.1101/2022.06.24.497523) (2022)


Documentation
* verkko on [GitHub](https://github.com/marbl/verkko)


Important Notes
* Module Name: verkko (see [the modules page](/apps/modules.html) for more information)
* verkko uses snakemake to run a complex workflow. For small data sets it can be run in local mode as a single job
 but for realistic size data it is likely to be used in cluster mode with `--slurm`
* Example files in `$VERKKO_TEST_DATA`



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run verkko in local mode for the smaller test
data. Don't move data to lscratch for slurm mode. Sample session:



```

[user@biowulf]$ **sinteractive --cpus-per-task=8 --mem=24g --gres=lscratch:100**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144]$ **module load verkko**
[user@cn3144]$ **cd /lscratch/$SLURM\_JOB\_ID**
[user@cn3144]$ **cp ${VERKKO\_TEST\_DATA:-none}/\* .**
[user@cn3144]$ **ls -lh**
total 363M
-rw-r--r-- 1 user group 119M Oct 26 12:07 hifi.fastq.gz
-rw-r--r-- 1 user group 245M Oct 26 12:07 ont.fastq.gz

[user@cn3144]$ **verkko -d asm --cleanup --local --local-memory 22 \
 --local-cpus $SLURM\_CPUS\_PER\_TASK --hifi ./hifi.fastq.gz --nano ./ont.fastq.gz**
Launching verkko branch  commit
Using snakemake 7.19.1
Building DAG of jobs...
Using shell: /usr/bin/bash
Provided cores: 8
Rules claiming more threads will be scaled down.
Provided resources: mem_gb=22
Job stats:
job                            count    min threads    max threads
---------------------------  -------  -------------  -------------
buildGraph                         1              4              4
[...snip...]
[user@cn3144]$ **ls -lh asm**
total 36M
drwxr-xr-x 2 user group 4.0K Oct 26 12:16 1-buildGraph
drwxr-xr-x 2 user group 4.0K Oct 26 12:26 2-processGraph
drwxr-xr-x 2 user group 4.0K Oct 26 12:26 3-align
drwxr-xr-x 2 user group 4.0K Oct 26 12:26 4-processONT
drwxr-xr-x 2 user group 4.0K Oct 26 12:26 5-untip
drwxr-xr-x 2 user group 4.0K Oct 26 12:18 6-layoutContigs
drwxr-xr-x 2 user group 4.0K Oct 26 12:26 7-consensus
-rw-r--r-- 1 user group 4.5M Oct 26 12:26 assembly.fasta
-rw-r--r-- 1 user group   40 Oct 26 12:26 assembly.hifi-coverage.csv
-rw-r--r-- 1 user group 3.3M Oct 26 12:26 assembly.homopolymer-compressed.gfa
-rw-r--r-- 1 user group 419K Oct 26 12:26 assembly.homopolymer-compressed.layout
-rw-r--r-- 1 user group  105 Oct 26 12:26 assembly.homopolymer-compressed.noseq.gfa
-rw-r--r-- 1 user group   41 Oct 26 12:26 assembly.ont-coverage.csv
-rw-r--r-- 1 user group    0 Oct 26 12:11 emptyfile
-rw-r--r-- 1 user group  28M Oct 26 12:16 hifi-corrected.fasta.gz
-rwxr-xr-x 1 user group  340 Oct 26 12:11 snakemake.sh
-rw-r--r-- 1 user group 2.9K Oct 26 12:11 verkko.yml

[user@cn3144]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
In this example verkko is run in slurm mode so the main job doesn't need many resources



```

#!/bin/bash
module load verkko/1.1
cp -r "${VERKKO_TEST_DATA:-none}" input
verkko -d asm --cleanup --local --local-memory 10 --slurm \
    --local-cpus $SLURM_CPUS_PER_TASK --hifi input/hifi.fastq.gz --nano input/ont.fastq.gz

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=2 --mem=10g verkko.sh
```







