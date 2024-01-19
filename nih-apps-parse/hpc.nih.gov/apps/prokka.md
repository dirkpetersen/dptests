

document.querySelector('title').textContent = 'prokka on Biowulf';
prokka on Biowulf


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


Prokka is a pipeline for rapidly annotating prokaryotic genomes. It produces
GFF3, GBK and SQN files that are ready for editing in Sequin and ultimately
submitted to Genbank/DDJB/ENA.


### References:


* T. Seeman. *Prokka: rapid prokaryotic genome annotation.*
 Bioinformatics 2014, 30:2068-2069.
 [Pubmed](https://www.ncbi.nlm.nih.gov/pubmed/24642063)  | 
 PMC  | 
 [Journal](http://bioinformatics.oxfordjournals.org/content/30/14/2068.long)


Documentation
* [Home page](http://www.vicbioinformatics.com/software.prokka.shtml)
* [Manual](https://github.com/tseemann/prokka/blob/master/README.md)
* [Tutorial](https://github.com/microgenomics/tutorials/blob/master/pangenome.md)


Important Notes
* Module Name: prokka (see [the modules page](/apps/modules.html) for more information)
* prokka is a multithreaded application



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=2g --cpus-per-task=4**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$
node$ **module load prokka**
node$ **prokka --listdb**
[16:11:09] Looking for databases in: /opt/anaconda/bin/../db
[16:11:09] * Kingdoms: Archaea Bacteria Mitochondria Viruses
[16:11:09] * Genera: Enterococcus Staphylococcus
[16:11:09] * HMMs: HAMAP
[16:11:09] * CMs: Bacteria Viruses
node$ **cp /usr/local/apps/prokka/TEST\_DATA/GCA\_000021185.1\_ASM2118v1\_genomic.fna .**
node$ **prokka --cpus 4 --force \
 --kingdom Bacteria \
 --outdir prokka\_GCA\_000021185 \
 --genus Listeria \
 --locustag GCA\_000021185 GCA\_000021185.1\_ASM2118v1\_genomic.fna**
[...snip...]

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. prokka.sh), which uses the input file 'prokka.in'. For example:



```

#! /bin/bash

function die {
    echo "$@" >&2
    exit 1
}

module load prokka/1.13 || die "Could not load prokka module"
cp /usr/local/apps/prokka/TEST_DATA/GCA_000021185.1_ASM2118v1_genomic.fna . \
    || die "Could not find test data"

prokka --cpus ${SLURM_CPUS_PER_TASK} --force \
    --kingdom Bacteria \
    --outdir prokka_GCA_000021185 \
    --genus Listeria \
    --locustag GCA_000021185 GCA_000021185.1_ASM2118v1_genomic.fna

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=6 --mem=3g --time=10 prokka.sh
```

This should create the following output directory:



```

./prokka_GCA_000021185
|-- GCA_000021185_10252016.err
|-- GCA_000021185_10252016.faa
|-- GCA_000021185_10252016.ffn
|-- GCA_000021185_10252016.fna
|-- GCA_000021185_10252016.fsa
|-- GCA_000021185_10252016.gbk
|-- GCA_000021185_10252016.gff
|-- GCA_000021185_10252016.log
|-- GCA_000021185_10252016.sqn
|-- GCA_000021185_10252016.tbl
`-- GCA_000021185_10252016.txt

```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Copy the example data



```

biowulf$ **cp -r /usr/local/apps/prokka/TEST\_DATA .**

```

Create a swarmfile (e.g. prokka.swarm). For example:



```

prokka --cpus ${SLURM_CPUS_PER_TASK} --force --kingdom Bacteria --outdir prokka_GCA_000008285 \
    --genus Listeria --locustag GCA_000008285 TEST_DATA/GCA_000008285.1_ASM828v1_genomic.fna
prokka --cpus ${SLURM_CPUS_PER_TASK} --force --kingdom Bacteria --outdir prokka_GCA_000021185 \
    --genus Listeria --locustag GCA_000021185 TEST_DATA/GCA_000021185.1_ASM2118v1_genomic.fna
prokka --cpus ${SLURM_CPUS_PER_TASK} --force --kingdom Bacteria --outdir prokka_GCA_000026705 \
    --genus Listeria --locustag GCA_000026705 TEST_DATA/GCA_000026705.1_ASM2670v1_genomic.fna
prokka --cpus ${SLURM_CPUS_PER_TASK} --force --kingdom Bacteria --outdir prokka_GCA_000168635 \
    --genus Listeria --locustag GCA_000168635 TEST_DATA/GCA_000168635.2_ASM16863v2_genomic.fna
prokka --cpus ${SLURM_CPUS_PER_TASK} --force --kingdom Bacteria --outdir prokka_GCA_000168815 \
    --genus Listeria --locustag GCA_000168815 TEST_DATA/GCA_000168815.1_ASM16881v1_genomic.fna
prokka --cpus ${SLURM_CPUS_PER_TASK} --force --kingdom Bacteria --outdir prokka_GCA_000196035 \
    --genus Listeria --locustag GCA_000196035 TEST_DATA/GCA_000196035.1_ASM19603v1_genomic.fna

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f prokka.swarm -g 2 -t 6 --module prokka/1.13
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g #  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t #  Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module prokka  Loads the prokka module for each subjob in the swarm 
 | |
 | |
 | |








