

document.querySelector('title').textContent = 'HMMRATAC on Biowulf';
HMMRATAC on Biowulf


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



A Hidden Markov ModeleR for ATAC-seq. HMMRATAC splits a single ATAC-seq dataset into nucleosome-free and nucleosome-enriched signals, learns the unique chromatin structure around accessible regions, and then predicts accessible regions across the entire genome.



### Reference:


* [Tarbell, Evan D., and Tao Liu. "HMMRATAC: a Hidden Markov ModeleR for ATAC-seq." *BioRxiv* (2019): 306621.](https://academic.oup.com/nar/article/47/16/e91/5519166)


Documentation
* [HMMRATAC GitHub repository](https://github.com/LiuLabUB/HMMRATAC)
* [HMMRATAC Guide](https://github.com/LiuLabUB/HMMRATAC/blob/master/HMMRATAC_Guide.md)


Important Notes
* Module Name: HMMRATAC (see [the modules page](/apps/modules.html) for more information)
 * Tests suggest that HMMRATAC will automatically try to use the number of threads that you allocate, but they are only utilized in short bursts. 
 * In testing, HMMRATAC used around 3GB of memory for each CPU allocated.
 * Environment variables set 
	+ HMMRATAC\_HOME* HMMRATAC is a java jar file and should be called like so:

	+ java -jar $HMMRATAC\_HOME/HMMRATAC.jar [options and arguments]



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):


First, use samtools to prepare a sample.



```

[user@biowulf ~]$ **sinteractive --cpus-per-task=4 --mem=12G --gres=lscratch:10**
salloc.exe: Pending job allocation 41216741
salloc.exe: job 41216741 queued and waiting for resources
salloc.exe: job 41216741 has been allocated resources
salloc.exe: Granted job allocation 41216741
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3125 are ready for job
srun: error: x11: no local DISPLAY defined, skipping

[user@cn3125 ~]$ **cd /data/$USER/HMMRATAC/testdata**

[user@cn3125 testdata]$ **ls**
ATACseq.bam

[user@cn3125 testdata]$ **module load samtools**
[+] Loading samtools 1.9  ...

[user@cn3125 testdata]$ **samtools sort -@4 -m 1800M -T /lscratch/$SLURM\_JOB\_ID/ATACseq.bam \
 -o ATACseq.sorted.bam ATACseq.bam**
[bam_sort_core] merging from 0 files and 4 in-memory blocks...

[user@cn3125 testdata]$ **samtools index -@4 ATACseq.sorted.bam ATACseq.sorted.bam.bai**

[user@cn3125 testdata]$ **samtools view -H ATACseq.sorted.bam | \
 perl -ne 'if(/^@SQ.\*?SN:(\w+)\s+LN:(\d+)/){print $1,"\t",$2,"\n"}' > \
 genome.info**

```

Then use HMMRATAC to analyze the sample. (Note that these --upper and --lower option/argument pairs may not be realistic.)



```

[user@cn3125 testdata]$ **module load HMMRATAC**
[+] Loading HMMRATAC  1.2.9  on cn3125
[+] Loading java 12.0.1  ...

[user@cn3125 testdata]$ **HMMRATAC\_HOME/HMMRATAC.jar --upper 100 --lower 2 \
 --bam ATACseq.sorted.bam \
 --index ATACseq.sorted.bam.bai \
 --genome genome.info**

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. HMMRATAC.sh). For example:



```

#!/bin/bash
#SBATCH --job-name="HMMRATAC-test"
#SBATCH --mem=12g
#SBATCH --ntasks=1
#SBATCH --partition=quick
#SBATCH --time=2:0:0
#SBATCH --cpus-per-task=4
#SBATCH --error=/data/$USER/HMMRATAC-test.e
#SBATCH --output=/data/$USER/HMMRATAC-test.o

module load HMMRATAC
cd /data/$USER/HMMRATAC/testdata

java -jar $HMMRATAC_HOME/HMMRATAC.jar \
    -b ATACseq.sorted.bam \
    -i ATACseq.sorted.bam.bai \
    -g genome.info

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
[user@biowulf ~]$ **sbatch HMMRATAC.sh**
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. HMMRATAC.swarm). For example:



```

java -jar $HMMRATAC_HOME/HMMRATAC.jar -b 1.sorted.bam -i 1.sorted.bam.bai -g 1.info
java -jar $HMMRATAC_HOME/HMMRATAC.jar -b 2.sorted.bam -i 2.sorted.bam.bai -g 2.info
java -jar $HMMRATAC_HOME/HMMRATAC.jar -b 3.sorted.bam -i 3.sorted.bam.bai -g 3.info
java -jar $HMMRATAC_HOME/HMMRATAC.jar -b 4.sorted.bam -i 4.sorted.bam.bai -g 4.info

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
[user@biowulf ~]$ **swarm -f HMMRATAC.swarm [-g #] [-t #] --module HMMRATAC**
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module HMMRATAC Loads the HMMRATAC module for each subjob in the swarm 
 | |
 | |
 | |








