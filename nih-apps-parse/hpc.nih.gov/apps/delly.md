

document.querySelector('title').textContent = 'delly on Biowulf';
delly on Biowulf


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


delly is an integrated structural variant prediction method that can detect
deletions, tandem duplications, inversions and translocations at
single-nucleotide resolution in short-read massively parallel sequencing data.
It uses paired-ends and split-reads to sensitively and accurately delineate
genomic rearrangements throughout the genome.


### References:


* Tobias Rausch, Thomas Zichner, Andreas Schlattl, Adrian M. Stuetz, 
 Vladimir Benes, Jan O. Korbel. *Delly: structural variant discovery by 
 integrated paired-end and split-read analysis.*. Bioinformatics 2012, 
 28:333-339.
 [Pubmed](https://www.ncbi.nlm.nih.gov/pubmed/22962449)  | 
 [PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3436805/)  | 
 [Journal](http://bioinformatics.oxfordjournals.org/content/28/18/i333.abstract)


Documentation
* [GitHub](https://github.com/tobiasrausch/delly)


Important Notes
* Module Name: delly (see [the modules page](/apps/modules.html) for more information)
* delly is a multithreaded application. The delly module sets the `$OMP_NUM_THREADS`
 environment variable automatically to match `$SLURM_CPUS_PER_TASK`. However,
 note that delly primarily parallelizes on the sample level, so there is no benefit to allocating
 multiple CPUs when processing a single sample.
* Example files in `$DELLY_TEST_DATA`
* Files with regions to exclude from calling for some genomes can be found in
 `$DELLY_EXCL_FILES`



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=5g**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load delly/0.7.8**
[user@cn3144 ~]$ **delly**
**********************************************************************
Program: Delly
This is free software, and you are welcome to redistribute it under
certain conditions (GPL); for license details use '-l'.
This program comes with ABSOLUTELY NO WARRANTY; for details use '-w'.

Delly (Version: 0.7.8)
Contact: Tobias Rausch (rausch@embl.de)
**********************************************************************

Usage: delly  

Commands:

 call discover and genotype structural variants
 merge merge structural variants across VCF/BCF files and within a single VCF/BCF file
 filter filter somatic or germline structural variants

[user@cn3144 ~]$ **cp $DELLY\_TEST\_DATA/\* .**
[user@cn3144 ~]$ # calling somatic SVs
[user@cn3144 ~]$ **delly call -o test.bcf -g DEL.fa DEL.bam**
[...snip...]
[user@cn3144 ~]$ **module load samtools**
[user@cn3144 ~]$ **bcftools view test.bcf**
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##fileDate=20180308
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=INV,Description="Inversion">
##ALT=<ID=BND,Description="Translocation">
##ALT=<ID=INS,Description="Insertion">
...

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. delly.sh) similar to the following example:



```

#! /bin/bash

function die {
    echo "$@" >&2
    exit 1
}

module load delly/0.7.8 || die "Could not load module"
cd /data/$USER/data_for_delly
delly call -o delly_calls.bcf -g ref.fa \
    sample1.bam sample2.bam sample3.bam sample4.bam

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=4 --mem=10g delly.sh
```

Loading the module as part of the batch script will automatically
set the OMP\_NUM\_THREADS variable to match the number of allocated
CPUs. If not loading the module in the batch script, please set
OMP\_NUM\_THREADS explicitly.


Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. delly.swarm). For example:



```

export OMP_NUM_THREADS=2; cd /data/$USER/dir1; delly call -o del.bcf -g ref.fa sample1.bam
export OMP_NUM_THREADS=2; cd /data/$USER/dir2; delly call -o del.bcf -g ref.fa sample2.bam
export OMP_NUM_THREADS=2; cd /data/$USER/dir3; delly call -o del.bcf -g ref.fa sample3.bam
export OMP_NUM_THREADS=2; cd /data/$USER/dir4; delly call -o del.bcf -g ref.fa sample4.bam

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f delly.swarm -g 10 -t 2 --module delly/0.7.8
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g #  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t #  Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module delly  Loads the delly module for each subjob in the swarm 
 | |
 | |
 | |








