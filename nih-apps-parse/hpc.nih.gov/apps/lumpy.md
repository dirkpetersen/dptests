

document.querySelector('title').textContent = 'Lumpy on Biowulf';
Lumpy on Biowulf


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



Lumpy detects structural variants from high throughput sequencing data
using read-pair, split read, and other signals in parallel.



There are two ways of using lumpy:



* `lumpyexpress`: Automated and standardized analysis for
 common situations
* `lumpy`: Flexible breakpoint detection for more
 complicated use cases


Either way, lumpy expects input alignments as created by `bwa
mem`.


### References:


* Ryan M. Layer, Colby Chiang, Aaron R. Quinlan and Ira M. Hall. 
 *LUMPY: a probabilistic framework for structural variant discovery*
 Genome Biol. 2014, 15:R84.
 [Pubmed](http://www.ncbi.nlm.nih.gov/pubmed/24970577) | 
 [PMC](http://www.ncbi.nlm.nih.gov/pmc/articles/pmid/24970577/)  | 
 [Journal](http://www.genomebiology.com/2014/15/6/R84)


Documentation
* [lumpy Main Site](https://github.com/arq5x/lumpy-sv)


Important Notes
* Module Name: lumpy (see [the modules page](/apps/modules.html) for more information)
* Multithreaded app
* environment variables set 
	+ $LUMPY\_HOME
	+ $LUMPY\_CONFIG



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

[user@cn3144 ~]$ **module load lumpy**

[user@cn3144 ~]$ **lumpyexpress -t2 -B /path/to/input.bam**
Sourcing executables from
/usr/local/apps/lumpy/0.2.11/bin/lumpyexpress.config ...

Checking for required python modules
(/usr/local/Anaconda/envs/py2.7.9/bin/python)...
samblaster: Version 0.1.22
samblaster: Inputting from stdin
samblaster: Outputting to stdout
samblaster: Opening gcat_set_053.bam.vcf.Kzd6atukx0F1/disc_pipe for write.
samblaster: Opening gcat_set_053.bam.vcf.Kzd6atukx0F1/spl_pipe for write.
samblaster: Loaded 25 header sequence entries.
samblaster: Output 27571 discordant read pairs to gcat_set_053.bam.vcf.Kzd6atukx0F1/disc_pipe
samblaster: Output 0 split reads to gcat_set_053.bam.vcf.Kzd6atukx0F1/spl_pipe
samblaster: Marked 0 of 39432108 (0.00%) read ids as duplicates using 1188k memory in 1M36S(96.300S) C
PU seconds and 19M47S(1187S) wall time.
Removed 5 outliers with isize >= 404
Running LUMPY... 
chrM    1000000
[...snip...]

[user@cn3144 ~]$ **lumpyexress -h**
usage:   lumpyexpress [options]

options:
     -B FILE  full BAM file(s) (comma separated) (required)
     -S FILE  split reads BAM file(s) (comma separated)
     -D FILE  discordant reads BAM files(s) (comma separated)
     -o FILE  output file [fullBam.bam.vcf]
     -x FILE  BED file to exclude
     -P       output probability curves for each variant
     -m INT   minimum sample weight for a call [4]
     -r FLOAT trim threshold [0]
     -T DIR   temp directory [./output_prefix.XXXXXXXXXXXX]
     -t N     number of threads [1]
     -k       keep temporary files

     -K FILE  path to lumpyexpress.config file
                (default: same directory as lumpyexpress)
     -v       verbose
     -h       show this message

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. lumpy.sh). For example:



```

#!/bin/bash

module load lumpy || exit 1
lumpyexpress -t $SLURM_CPUS_PER_TASK -B input.bam -o output.vcf

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=10 --mem=10G lumpy.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. lumpy.swarm). For example:



```

lumpyexpress -t $SLURM_CPUS_PER_TASK -B input1.bam -o output1.vcf
lumpyexpress -t $SLURM_CPUS_PER_TASK -B input2.bam -o output2.vcf
lumpyexpress -t $SLURM_CPUS_PER_TASK -B input3.bam -o output3.vcf

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f lumpy.swarm -g 10 -t 10 --module lumpy
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module lumpy Loads the lumpy module for each subjob in the swarm 
 | |
 | |
 | |








