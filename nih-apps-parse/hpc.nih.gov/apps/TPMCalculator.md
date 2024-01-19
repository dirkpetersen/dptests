

document.querySelector('title').textContent = 'TPMCalculator on Biowulf';
TPMCalculator on Biowulf


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



TPMCalculator quantifies mRNA abundance directly from the alignments by parsing BAM files. The input parameters are the same GTF files used to generate the alignments, and one or multiple input BAM file(s) containing either single-end or paired-end sequencing reads. The TPMCalculator output is comprised of four files per sample reporting the TPM values and raw read counts for genes, transcripts, exons and introns respectively.



### References:


* Roberto Vera Alvarez, Lorinc Sandor Pongor, Leonardo Mariño-Ramírez, David Landsman; TPMCalculator: one-step software to quantify mRNA abundance of genomic features, Bioinformatics, , bty896, <https://doi.org/10.1093/bioinformatics/bty896>


Documentation
* [TPMCalculator on github](https://github.com/ncbi/TPMCalculator)
* [TPMCalculator Wiki](https://github.com/ncbi/TPMCalculator/wiki)


Important Notes
* Module Name: TPMCalculator (see [the modules page](/apps/modules.html) for more information)
 * Singlethreaded.
* Example files in /usr/local/apps/TPMCalculator/TEST\_DATA* Reference data in fdb/igenomes/



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

[user@cn3144 ~]$ **module load TPMCalculator**

[user@cn3144 ~]$ **cp -r /usr/local/apps/TPMCalculator/TEST\_DATA /data/$USER**

[user@cn3144 ~]$ **cd /data/$USER/TEST\_DATA**

[user@cn3144 ~]$ **TPMCalculator -v -g /fdb/igenomes/Homo\_sapiens/UCSC/hg38/Annotation/Genes/genes.gtf -b SRR2126790\_sorted.bam -q 255**
Reading GTF file ...
Done in 57.6513 seconds
Parsing sample: SRR2126790_sorted
 3087949 reads processed in 80.0002 seconds
Printing results
Total time: 139.753 seconds

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. TPMCalculator.sh). For example, the following script will copy the sample data, run the program, and compare the output to 
the validated results. 



```

#!/bin/bash
set -e
module load TPMCalculator

cp -r /usr/local/apps/TPMCalculator/TEST_DATA  /data/$USER

cd /data/$USER/TEST_DATA

./run.sh

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] TPMCalculator.sh
```

Note: the test job in the above batch script uses less than 1 GB of memory, so the default allocated memory is sufficient. For your own jobs, if they need more memory, you may need to use the --mem flag.

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. TPMCalculator.swarm). For example:



```

TPMCalculator -v -g /fdb/igenomes/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.gtf   -b file1.bam   -q 255
TPMCalculator -v -g /fdb/igenomes/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.gtf   -b file2.bam   -q 255
TPMCalculator -v -g /fdb/igenomes/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.gtf   -b file3.bam   -q 255
[...]

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f TPMCalculator.swarm [-g #] --module TPMCalculator
```

where


|  |  |  |  |
| --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | --module TPMCalculator Loads the TPMCalculator module for each subjob in the swarm 
 | |
 | |








