

document.querySelector('title').textContent = 'Sequenza-utils on Biowulf';
Sequenza-utils on Biowulf


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



Sequenza-utils is The supporting python library for the sequenza R package.



Documentation
* [Sequenza-utils Main Site](https://sequenza-utils.readthedocs.io/en/latest/)


Important Notes
* Module Name: sequenza-utils (see [the modules page](/apps/modules.html) for more information)
* singlethreaded app



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load sequenza-utils**

[user@cn3144 ~]$ **sequenza-utils -h**
usage: sequenza-utils [-h] [-v]
                      {bam2seqz,gc_wiggle,pileup2acgt,seqz_binning,snp2seqz}
                      ...

Sequenza Utils is an ensemble of tools capable of perform various tasks, primarily aimed to convert bam/pileup files to a format usable by the sequenza R package

positional arguments:
    bam2seqz        Process a paired set of BAM/pileup files (tumor and
                    matching normal), and GC-content genome-wide
                    information, to extract the common positions withA and
                    B alleles frequencies
    gc_wiggle       Given a fasta file and a window size it computes the GC
                    percentage across the sequences, and returns a file in
                    the UCSC wiggle format.
    pileup2acgt     Parse the format from the samtools mpileup command, and
                    report the occurrence of the 4 nucleotides in each
                    position.
    seqz_binning    Perform the binning of the seqz file to reduce file
                    sizeand memory requirement for the analysis.
    snp2seqz        Parse VCFs and other variant and coverage formats to
                    produce seqz files

optional arguments:
  -h, --help        show this help message and exit
  -v, --verbose     Show all logging information

[user@cn3144 ~]$ **sequenza-utils gc\_wiggle -f /fdb/app\_testdata/fasta/R64-1-1.cdna\_nc.fa -o sequenza.out**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. sequenza-utils.sh). For example:



```

#!/bin/bash
module load sequenza-utils
sequenza-utils gc_wiggle -f /fdb/app_testdata/fasta/R64-1-1.cdna_nc.fa -o sequenza.out
```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch sequenza-utils.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. sequenza-utils.swarm). For example:



```

sequenza-utils gc_wiggle -f s1.fa -o s1.out
sequenza-utils gc_wiggle -f s2.fa -o s2.out
sequenza-utils gc_wiggle -f s3.fa -o s3.out

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f sequenza-utils.swarm --module sequenza-utils
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module sequenza-utils Loads the sequenza-utils module for each subjob in the swarm 
 | |
 | |
 | |








