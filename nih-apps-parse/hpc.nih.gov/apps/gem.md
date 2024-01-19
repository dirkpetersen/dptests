

document.querySelector('title').textContent = 'Gem on HPC';
Gem on HPC


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

  GEM: High resolution peak calling and motif discovery for ChIP-seq and 
 ChIP-exo data.


GEM is a scientific software for studying protein-DNA interaction at high 
 resolution using ChIP-seq/ChIP-exo data. It can also be applied to CLIP-seq 
 and Branch-seq data.   

 GEM links binding event discovery and motif discovery with positional priors 
 in the context of a generative probabilistic model of ChIP data and genome 
 sequence, resolves ChIP data into explanatory motifs and binding events 
 at unsurpassed spatial resolution. GEM reciprocally improves motif discovery 
 using binding event locations, and binding event predictions using discovered 
 motifs.


GEM has following features:


* Exceptionally high spatial resolution on binding event prediction (aka 
 peak calling)
* Highly accurate de novo motif discovery
* Resolves closely spaced (less than 500bp) homotypic events that appear 
 as a single cluster of reads
* Enables analysis of spatial binding constraints of multiple transcription 
 factors, for predicting TF dimer binding sites, enhanceosomes, etc.
* Analyzes ChIP-seq, ChIP-exo, CLIP-seq and Branch-seq data, single-end 
 or paired-end
* Runs in single-condition mode or multi-condition mode


### References:

 * <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3415389/>


Documentation * <http://groups.csail.mit.edu/cgs/gem/>



Important Notes * Module Name: gem (see [the modules 
 page](/apps/modules.html) for more information)
* Multithreaded
* environment variables set
	+ `$GEMJARPATH`
	+ `$GEMJAR`
	+ $GEM\_JAR
	+ $GEM\_JARPATH
* Example files in /usr/local/apps/gem/TEST\_DATA





Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive -c 10 --mem 10g**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144]$ **module load gem**
[user@cn3144]$ **cp -r ${GEM\_TEST\_DATA:-none}/ .
[user@cn3144]$ **java -Xmx10g -jar $GEMJAR --t 8 --d Read\_Distribution\_default.txt \
--g mm10.chrom.sizes \
--genome /fdb/igenomes/Mus\_musculus/UCSC/mm10/Sequence/Chromosomes/ \
--s 2000000000 --expt SRX000540\_mES\_CTCF.bed --ctrl SRX000543\_mES\_GFP.bed \
--f BED --out mouseCTCF --k\_min 6 --k\_max 13****
```


```

[user@cn3144]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```




Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. batch.sh). For example:



```

#!/bin/bash
set -e
module load gem
java -Xmx10g -jar $GEMJAR --t 8 --d $GEM_TEST_DATA/Read_Distribution_default.txt \  
--g $GEM_TEST_DATA/mm10.chrom.sizes \  
--genome /fdb/igenomes/Mus_musculus/UCSC/mm10/Sequence/Chromosomes/ \  
--s 2000000000 --expt SRX000540_mES_CTCF.bed --ctrl SRX000543_mES_GFP.bed \  
--f BED --out mouseCTCF --k_min 6 --k_max 13
```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=10 --mem=10g batch.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. job.swarm). For example:



```

cd dir1; gem command
cd dir2; gem command
cd dir3; gem command

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f job.swarm -g 10 -t 4 --module gem
```

where
 

|  |  |
| --- | --- |
| -g *#*  | Number of Gigabytes of memory required for each process (1 line in the swarm command file)
  |
| -t *#* | Number of threads/CPUs required for each process (1 line in the swarm command file).
  |
| --module  | Loads the module for each subjob in the swarm  |




