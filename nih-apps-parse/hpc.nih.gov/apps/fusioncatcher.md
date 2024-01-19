

document.querySelector('title').textContent = 'fusioncatcher on Biowulf';
fusioncatcher on Biowulf


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



FusionCatcher searches for novel/known somatic fusion genes, translocations, 
and chimeras in RNA-seq data (paired-end or single-end reads from Illumina NGS 
platforms like Solexa/HiSeq/NextSeq/MiSeq) from diseased samples.



### References:


* [Nicorici, Daniel, et al. "FusionCatcher-a tool for finding somatic fusion genes in paired-end RNA-sequencing data." *bioRxiv* (2014): 011650.](http://www.biorxiv.org/content/early/2014/11/19/011650)
* [Tuna, Musaffe. "Next‐Generation Sequencing in Cancer: Tools for Fusion Gene Detection." *eLS* (2015).](http://onlinelibrary.wiley.com/doi/10.1002/9780470015902.a0025848/full)


Documentation
* [GitHub](https://github.com/ndaniel/fusioncatcher)
* [sourceforge](http://sourceforge.net/projects/fusioncatcher/)
* [fusioncatcher manual](https://github.com/ndaniel/fusioncatcher/blob/master/doc/manual.md)


Important Notes
* Module Name: fusioncatcher (see [the modules page](/apps/modules.html) for more information)
* Reference data in $FUSIONCATCHER\_REFERENCE\_DATA
**Disk Usage**  

FusionCatcher jobs may use a great deal of disk space for temporary files. In
some instances disk usage may run into the hundreds of GB range. If the 
--output directory is located on network storage, the I/O load can lead to poor
performance and affect other users. Local scratch space (lscratch) 
should therefore be used when jobs are submitted to the batch system (as in the 
examples below). See [the Biowulf user guide](https://hpc.nih.gov/docs/userguide.html#local) for more information about using lscratch.



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive --mem 20g**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load fusioncatcher**
[user@cn3144 ~]$ **mkdir -p /data/$USER/fusioncatcher/test**
[user@cn3144 ~]$ **cd /data/$USER/fusioncatcher/test**
[user@cn3144 test]$ **wget http://sourceforge.net/projects/fusioncatcher/files/test/reads\_{1,2}.fq.gz**
[user@cn3144 ~]$ **cd ..**
[user@cn3144 ~]$ **fusioncatcher \
 --input ./test/ \
 --output ./test-results/ \
 --threads=2**
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. fusioncatcher.sh). For example:



```

#!/bin/sh
set -e

module load fusioncatcher
fusioncatcher \
 --input /data/$USER/fusioncatcher/test/ \
 --output /lscratch/$SLURM_JOB_ID/test-results/ \
 --threads=2

mkdir /data/$USER/fusioncatcher/$SLURM_JOB_ID
mv /lscratch/$SLURM_JOB_ID/test-results /data/$USER/fusioncatcher/$SLURM_JOB_ID/test-results

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] --mem 20g --gres lscratch:10 --time 30 fusioncatcher.sh
```


Note that larger jobs may need more memory, lscratch space, and longer walltime.
A typical single-sample paired end RNA-Seq run should requred less than 50GB memory and 100GB lscratch space.



Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.

FusionCatcher has a built-in "batch mode" which can be accessed using the 
fusioncatcher-batch command. This command accepts a filename as input. The
file should be tab delimited text listing paired input and output. Interested 
users should see
[the FusionCatcher manual](https://github.com/ndaniel/fusioncatcher/blob/master/doc/manual.md)
for more details. However, users should note that the batch scripts used for
fusioncatcher-batch jobs can easily be converted to swarm scripts and run in 
parallel for great speed increases. 




In this section, we will begin by considering a fusioncatcher-batch job and see
how to convert it to a swarm job for increased efficiency. An example 
FusionCatcher batch script that uses ftp URLs as input can be downloaded 
[here](http://sourceforge.net/projects/fusioncatcher/files/examples/illumina-bodymap2.txt).
For brevity we will refer to the contents of this file in abbreviated form
below.




The batch script looks something like so:




```

ftp.uk/72.fastq   thyroid
ftp.uk/73.fastq   testis
ftp.uk/74.fastq   ovary
ftp.uk/75.fastq   leukocyte
ftp.uk/76.fastq   skeletal muscle
ftp.uk/77.fastq   prostate
ftp.uk/78.fastq   lymph node
ftp.uk/79.fastq   lung
ftp.uk/80.fastq   adipose
ftp.uk/81.fastq   adrenal
ftp.uk/82.fastq   brain
ftp.uk/83.fastq   breast
ftp.uk/84.fastq   colon
ftp.uk/85.fastq   kidney
ftp.uk/86.fastq   heart
ftp.uk/87.fastq   liver 

```


The ftp URLs on the left denote input for the batch script, and the names on 
the right give the locations of directories that should contain output when the
script is run. 




This syntax can be adapted to make a swarm file:




```

fusioncatcher -i ftp.uk/72.fastq -o /lscratch/$SLURM_JOB_ID/thyroid         -p 2; \
    mv /lscratch/$SLURM_JOB_ID/thyroid /data/$USER/thyroid 
fusioncatcher -i ftp.uk/73.fastq -o /lscratch/$SLURM_JOB_ID/testis          -p 2; \
    mv /lscratch/$SLURM_JOB_ID/thyroid /data/$USER/testis
fusioncatcher -i ftp.uk/74.fastq -o /lscratch/$SLURM_JOB_ID/ovary           -p 2; \
    mv /lscratch/$SLURM_JOB_ID/thyroid /data/$USER/ovary
fusioncatcher -i ftp.uk/75.fastq -o /lscratch/$SLURM_JOB_ID/leukocyte       -p 2; \
    mv /lscratch/$SLURM_JOB_ID/thyroid /data/$USER/leukocyte
fusioncatcher -i ftp.uk/76.fastq -o /lscratch/$SLURM_JOB_ID/skeletal_muscle -p 2; \
    mv /lscratch/$SLURM_JOB_ID/thyroid /data/$USER/skeletal_muscle
fusioncatcher -i ftp.uk/77.fastq -o /lscratch/$SLURM_JOB_ID/prostate        -p 2; \
    mv /lscratch/$SLURM_JOB_ID/thyroid /data/$USER/prostate
fusioncatcher -i ftp.uk/78.fastq -o /lscratch/$SLURM_JOB_ID/lymph_node      -p 2; \
    mv /lscratch/$SLURM_JOB_ID/thyroid /data/$USER/lymph_node
fusioncatcher -i ftp.uk/79.fastq -o /lscratch/$SLURM_JOB_ID/lung            -p 2; \
    mv /lscratch/$SLURM_JOB_ID/thyroid /data/$USER/lung
fusioncatcher -i ftp.uk/80.fastq -o /lscratch/$SLURM_JOB_ID/adipose         -p 2; \
    mv /lscratch/$SLURM_JOB_ID/thyroid /data/$USER/adipose
fusioncatcher -i ftp.uk/81.fastq -o /lscratch/$SLURM_JOB_ID/adrenal         -p 2; \
    mv /lscratch/$SLURM_JOB_ID/thyroid /data/$USER/adrenal
fusioncatcher -i ftp.uk/82.fastq -o /lscratch/$SLURM_JOB_ID/brain           -p 2; \
    mv /lscratch/$SLURM_JOB_ID/thyroid /data/$USER/brain
fusioncatcher -i ftp.uk/83.fastq -o /lscratch/$SLURM_JOB_ID/breast          -p 2; \
    mv /lscratch/$SLURM_JOB_ID/thyroid /data/$USER/breast
fusioncatcher -i ftp.uk/84.fastq -o /lscratch/$SLURM_JOB_ID/colon           -p 2; \
    mv /lscratch/$SLURM_JOB_ID/thyroid /data/$USER/colon
fusioncatcher -i ftp.uk/85.fastq -o /lscratch/$SLURM_JOB_ID/kidney          -p 2; \
    mv /lscratch/$SLURM_JOB_ID/thyroid /data/$USER/kidney
fusioncatcher -i ftp.uk/86.fastq -o /lscratch/$SLURM_JOB_ID/heart           -p 2; \
    mv /lscratch/$SLURM_JOB_ID/thyroid /data/$USER/heart
fusioncatcher -i ftp.uk/87.fastq -o /lscratch/$SLURM_JOB_ID/liver           -p 2; \
    mv /lscratch/$SLURM_JOB_ID/thyroid /data/$USER/liver

```


Note the substitution of underscores for spaces in the output directories.
Note also the -p argument limiting each sub-job to 2 processors. This is
appropriate because SLURM will assign each job 2 cpus (1 hyperthreaded core) by
default. Assuming that this file is saved as fusioncatcher.swarm, it could be 
submitted using the [swarm](/apps/swarm.html) command as follows:




```

[user@helix ~]$ **swarm -f fusioncatcher.swarm [-t #] -g 20 --time 20 --module fusioncatcher --gres=lscratch:10**

```


where


|  |  |  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module fusioncatcher Loads the fusioncatcher module for each subjob in the swarm 
 | --gres lscratch:10 Allocates 10 GB of local scratch space for each subjob in the swarm
 | |
 | |
 | |
 | |








