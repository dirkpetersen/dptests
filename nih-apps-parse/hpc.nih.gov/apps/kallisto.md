

document.querySelector('title').textContent = 'kallisto on Biowulf';
kallisto on Biowulf


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



kallisto is a program for quantifying abundances of transcripts from RNA-Seq data, or more generally of target sequences using high-throughput sequencing reads. It is based on the novel idea of pseudoalignment for rapidly determining the compatibility of reads with targets, without the need for alignment.



### References:


* NL Bray, H Pimentel, P Melsted and L Pachter, Near optimal probabilistic RNA-seq quantification, Nature Biotechnology 34, p 525--527 (2016). doi:[10.1038/nbt.3519](https://doi.org/10.1038/nbt.3519)


Documentation
* [Kallisto Main Site](http://pachterlab.github.io/kallisto/)
* [Google Group](https://groups.google.com/d/forum/kallisto-sleuth-users)


Important Notes
* Module Name: kallisto (see [the modules page](/apps/modules.html) for more information)
* Example files in /fdb/kallisto



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **cp -r /fdb/kallisto/\* .** 
[user@cn3144 ~]$ **kallisto index -i transcripts.kidx transcripts.fasta.gz**

[build] loading fasta file transcripts.fasta.gz
[build] k-mer length: 31
[build] counting k-mers ... done.
[build] building target de Bruijn graph ...  done 
[build] creating equivalence classes ...  done
[build] target de Bruijn graph has 27 contigs and contains 22118 k-mers 

[user@cn3144 ~]$ **kallisto quant \
 -i transcripts.kidx \
 -b 30 \
 -o kallisto\_out \
 --genomebam \
 --gtf transcripts.gtf.gz \
 --chromosomes chrom.txt \
 reads\_1.fastq.gz reads\_2.fastq.gz**
[quant] fragment length distribution will be estimated from the data
[index] k-mer length: 31
[index] number of targets: 14
[index] number of k-mers: 22,118
[index] number of equivalence classes: 20
Warning: 13 transcripts were defined in GTF file, but not in the index
[quant] running in paired-end mode
[quant] will process pair 1: reads_1.fastq.gz
                             reads_2.fastq.gz
[quant] finding pseudoalignments for the reads ... done
[quant] processed 10,000 reads, 9,413 reads pseudoaligned
[quant] estimated average fragment length: 178.02
[   em] quantifying the abundances ... done
[   em] the Expectation-Maximization algorithm ran for 52 rounds
[bstrp] running EM for the bootstrap: 30
[  bam] writing pseudoalignments to BAM format .. done
[  bam] sorting BAM files .. done
[  bam] indexing BAM file .. done

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. kallisto.sh). For example:



```

#!/bin/sh
set -e
module load kallisto

kallisto index -i transcripts.kidx transcripts.fasta.gz

kallisto quant \
 -i transcripts.kidx \
 -b 30 \
 -o kallisto_out \
 --genomebam \
 --gtf transcripts.gtf.gz \
 --chromosomes chrom.txt \
 reads_1.fastq.gz reads_2.fastq.gz

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] kallisto.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. kallisto.swarm). For example:



```

kallisto quant -i transcripts.kidx -b 30 -o kallisto/sample-01 --genomebam --gtf transcripts.gtf.gz --chromosomes chrom.txt sample01/reads_{1,2}.fastq.gz
kallisto quant -i transcripts.kidx -b 30 -o kallisto/sample-02 --genomebam --gtf transcripts.gtf.gz --chromosomes chrom.txt sample02/reads_{1,2}.fastq.gz
kallisto quant -i transcripts.kidx -b 30 -o kallisto/sample-03 --genomebam --gtf transcripts.gtf.gz --chromosomes chrom.txt sample03/reads_{1,2}.fastq.gz
kallisto quant -i transcripts.kidx -b 30 -o kallisto/sample-04 --genomebam --gtf transcripts.gtf.gz --chromosomes chrom.txt sample04/reads_{1,2}.fastq.gz

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f kallisto.swarm [-g #] [-t #] --module kallisto
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module kallisto Loads the kallisto module for each subjob in the swarm 
 | |
 | |
 | |








