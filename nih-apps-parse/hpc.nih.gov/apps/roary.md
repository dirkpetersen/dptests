

document.querySelector('title').textContent = 'roary on Biowulf';
roary on Biowulf


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


 Roary is a high speed pan genome pipeline, which takes annotated assemblies
in GFF3 format (produced by 
[Prokka](https://hpc.nih.gov/apps/prokka.html)) and calculates the pan
genome.


### References:


* Andrew J. Page et al., *Roary: Rapid large-scale prokaryote pan genome analysis*. 
 Bioinformatics 2015, 31:3691-3693.
 [Pubmed](https://www.ncbi.nlm.nih.gov/pubmed/26198102)  | 
 [PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4817141/)  | 
 [Journal](http://bioinformatics.oxfordjournals.org/content/31/22/3691.long)


Documentation
* [Home page](http://sanger-pathogens.github.io/Roary/)
* [GitHub](https://github.com/sanger-pathogens/Roary)
* [Tutorial](https://github.com/microgenomics/tutorials/blob/master/pangenome.md)


Important Notes
* Module Name: roary (see [the modules page](/apps/modules.html) for more information)
* roary is a multithreaded application. Make sure to match the number of cpus requested with the
 number of threads.
* Example files in `ROARY_TEST_DATA`



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
For this [interactive session](/docs/userguide.html#int) we will follow the
roary [Tutorial](https://github.com/microgenomics/tutorials/blob/master/pangenome.md). In addition
to roary we will also need prokka for annotating the bacterial genomes



```

[user@biowulf]$ **sinteractive --cpus-per-task=6 --mem=20g --gres=lscratch:20**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144]$ **module load roary prokka**
[user@cn3144]$ **cd /lscratch/$SLURM\_JOB\_ID**
[user@cn3144]$ **cp ${ROARY\_TEST\_DATA:-none}/\* .**
[user@cn3144]$ **ls -1**
GCA_000008285.1_ASM828v1_genomic.fna
GCA_000021185.1_ASM2118v1_genomic.fna
GCA_000026705.1_ASM2670v1_genomic.fna
GCA_000168635.2_ASM16863v2_genomic.fna
GCA_000168815.1_ASM16881v1_genomic.fna
GCA_000196035.1_ASM19603v1_genomic.fna

```

Annotate the genomes with prokka and create a pan genome with roary.



```

[user@cn3144]$ **fastas=( \*.fna )**
[user@cn3144]$ **for fasta in ${fastas[@]}; do
 acc=${fasta:0:13}
 prokka --kingdom Bacteria --outdir prokka\_$acc --genus Listeria \
 --locustag $acc --prefix $acc --cpus $SLURM\_CPUS\_PER\_TASK $fasta
 done**
[user@cn3144]$ **roary -p $SLURM\_CPUS\_PER\_TASK -f ./roary\_demo -e -n -v -r \*/\*.gff**
[...snip...]
[user@cn3144]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. roary.sh), which uses the input file 'roary.in'. For example:



```

#! /bin/bash
# this file is roary.batch

function die {
    echo "$@" >&2
    exit 1
}

wd=$PWD
module load roary/3.12.0 || die "Could not load modules"
cd /lscratch/$SLURM_JOB_ID || die "no lscratch"
roary -p ${SLURM_CPUS_PER_TASK} \
  -f roary_out -e -n -r -v ./gff/*.gff \
&& mv roary_out $wd

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --gres=lscratch:10 --cpus-per-task=6 --mem=6g roary.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. roary.swarm). For example:



```

roary -p $SLURM_CPUS_PER_TASK -f ./roary_species1 -e -n -v -r species1/*.gff
roary -p $SLURM_CPUS_PER_TASK -f ./roary_species2 -e -n -v -r species2/*.gff
roary -p $SLURM_CPUS_PER_TASK -f ./roary_species3 -e -n -v -r species3/*.gff

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f roary.swarm -g 6 -t 6 --module roary/3.12.0
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g #  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t #  Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module roary  Loads the roary module for each subjob in the swarm 
 | |
 | |
 | |








