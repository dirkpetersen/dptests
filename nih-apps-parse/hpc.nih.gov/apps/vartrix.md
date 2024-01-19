

document.querySelector('title').textContent = 'Vartrix on Biowulf';
Vartrix on Biowulf


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



VarTrix is a software tool for extracting single cell variant information from 10x Genomics single cell data. VarTrix will take a set of previously defined variant calls and use that to identify those variants in the single cell data.



### Web site


* [Home page](https://github.com/10XGenomics/vartrix)


Documentation
* [Vatrix Documentation](https://github.com/10XGenomics/vartrix)


Important Notes
* Module Name: vartrix (see [the modules page](/apps/modules.html) for more information)



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

[user@cn3144 ~]$ **module load vartrix**

[user@cn3144 ~]$ **cd /data/${USER}**

[user@cn3144 ~]$ **cp ${VARTRIX\_TEST\_DATA}/\* .**

[user@cn3144 ~]$  **vartrix --bam test\_dna.bam --cell-barcodes dna\_barcodes.tsv --fasta test\_dna.fa --out-matrix test\_dna.out --vcf test\_dna.vcf**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. vartrix\_job.sh). For example:



```

#!/bin/bash
set -e
module load vartrix
vartrix --bam test_dna.bam \
   --cell-barcodes dna_barcodes.tsv \
   --fasta test_dna.fa \
   --out-matrix test_dna.out \
   --vcf test_dna.vcf

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] vartrix_job.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. vartrix\_jobs.swarm). For example:



```

vartrix --bam test_dna1.bam \ 
   --cell-barcodes dna_barcodes1.tsv \
   --fasta test_dna1.fa \
   --out-matrix test_dna1.out \
   --vcf test_dna1.vcf
vartrix --bam test_dna2.bam \
   --cell-barcodes dna_barcodes2.tsv \
   --fasta test_dna2.fa \
   --out-matrix test_dna2.out \
   --vcf test_dna2.vcf
vartrix --bam test_dna3.bam \
   --cell-barcodes dna_barcodes3.tsv \
   --fasta test_dna3.fa \
   --out-matrix test_dna3.out \
   --vcf test_dna3.vcf

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f vartrix_jobs.swarm [-g #] [-t #] --module vartrix
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module vartrix Loads the vartrix module for each subjob in the swarm 
 | |
 | |
 | |








