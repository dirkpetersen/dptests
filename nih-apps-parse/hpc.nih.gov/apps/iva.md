

document.querySelector('title').textContent = 'iva on Biowulf';
iva on Biowulf


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


 Integrated virus assembler (IVA) is a de novo assembler designed for virus
genomes without repeat sequences, using Illumina read pairs from mixed
populations at extremely high and variable depth.


### References:


* Martin Hunt *et al.*. *IVA: accurate de novo assembly of 
 RNA virus genomes*. Bioinformatics 2015(31): 2374-2376.
 [Pubmed](http://www.ncbi.nlm.nih.gov/pubmed/25725497) | 
 [PMC](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4495290/) | 
 [Journal](http://bioinformatics.oxfordjournals.org/content/31/14/2374)


Documentation
* [Home page](http://sanger-pathogens.github.io/iva/)
* [Wiki](https://github.com/sanger-pathogens/iva/wiki)
* [GitHub](https://github.com/sanger-pathogens/iva)


Important Notes
* Module Name: iva (see [the modules page](/apps/modules.html) for more information)
* iva is a multithreaded application. Please match the number of allocated CPUs to the
 number of threads
* some iva versions mis-report their version
* some samtools errors may be spurious



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the built in test data
set



```

[user@biowulf]$ **sinteractive --cpus-per-task=2 --mem=6g**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144]$ **module load iva**
[user@cn3144]$ **iva --test --threads=2 --trimmomatic=$TRIMMOJAR test**
Running iva in test mode...
Copied input test files into here: PREFIX/test
Current working directory: PREFIX/test
Running iva on the test data with the command:
/usr/local/lib/python3.6/dist-packages/iva-1.0.9-py3.6.egg/EGG-INFO/scripts/iva --threads 2 \
    --trimmomatic /Trimmomatic-0.38/trimmomatic-0.38.jar --pcr_primers hiv_pcr_primers.fa -f reads_1.fq.gz -r reads_2.fq.gz iva.out
Finished running iva
Looks OK. Final output contigs file is: PREFIX/test/iva.out/contigs.fasta
[user@cn3144]$ **tree test**
test
|-- [user    672]  hiv_pcr_primers.fa
|-- [user   9.7K]  iva_contigs_no_trimmomatic.fasta
|-- [user   9.0K]  iva_contigs_with_trimmomatic.fasta
|-- [user   4.0K]  iva.out
|   |-- [user   3.1K]  adapters.fasta
|   |-- [user   9.8K]  contigs.fasta
|   `-- [user    369]  info.txt
|-- [user   3.6M]  reads_1.fq.gz
|-- [user   4.4M]  reads_2.fq.gz
`-- [user   9.0K]  reference.fasta

```

Run IVA with trimmomatic on the data copied to test by the automated 
test



```

[user@cn3144]$ **cd test**
[user@cn3144]$ **iva --threads 2 --pcr\_primers hiv\_pcr\_primers.fa \
 -f reads\_1.fq.gz -r reads\_2.fq.gz \
 --trimmomatic $TRIMMOJAR \
 iva.out2**  

[user@cn3144]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. iva.sh), which uses the input file 'iva.in'. For example:



```

#! /bin/bash

function fail {
  echo >2 "$@"
  exit 1
}

module load iva/1.0.11 || fail "could not load iva module"
iva --threads $SLURM_CPUS_PER_TASK \
  -f read1.fq.gz -r read2.fq.gz \
  --trimmomatic=$TRIMMOJAR \
  sample.iva || fail "iva return non-zero exit status"

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=6 iva.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. iva.swarm). For example:



```

iva -r sample1_r1.fq.gz -f sample1_r2.fq.gz --trimmomatic=$TRIMMOJAR --threads $SLURM_CPUS_PER_TASK sample1_out
iva -r sample2_r1.fq.gz -f sample2_r2.fq.gz --trimmomatic=$TRIMMOJAR --threads $SLURM_CPUS_PER_TASK sample2_out
iva -r sample3_r1.fq.gz -f sample3_r2.fq.gz --trimmomatic=$TRIMMOJAR --threads $SLURM_CPUS_PER_TASK sample3_out

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f iva.swarm -g 10 -t 4 --module iva/1.0.3
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g #  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t #  Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module iva  Loads the iva module for each subjob in the swarm 
 | |
 | |
 | |








