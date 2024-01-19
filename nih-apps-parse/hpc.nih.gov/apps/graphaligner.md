

document.querySelector('title').textContent = 'graphaligner on Biowulf';
graphaligner on Biowulf


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



Seed-and-extend program for aligning long error-prone reads to genome graphs.



### References:


* M. Rautiainen, T. Marschall. *GraphAligner: rapid and versatile sequence-to-graph alignment*.
 Genome Biology (2020). [PubMed](https://pubmed.ncbi.nlm.nih.gov/32972461/) | 
 [PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/pmid/32972461/) | 
 [Journal](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02157-2)


Documentation
* graphaligner on [GitHub](https://github.com/maickrau/GraphAligner)


Important Notes
* Module Name: graphaligner (see [the modules page](/apps/modules.html) for more information)
* GraphAligner is a multithreaded application. Please match the number of threads with the number of allocated
 CPUs
* Example files in `$GRAPHALIGNER_TEST_DATA`



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --cpus-per-task=2 --gres=lscratch:10**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144]$ **module load graphaligner**
[user@cn3144]$ **cp -r ${GRAPHALIGNER\_TEST\_DATA:-none} test**
[user@cn3144]$ **GraphAligner -t $SLURM\_CPUS\_PER\_TASK -g test/graph.gfa -f test/read.fa -a test/aln.gaf -x vg**
[user@cn3144]$ **cat test/aln.gaf**
read    71      0       71      +       >1>2>4  87      3       73      67      72      60      NM:i:5  AS:f:56.3       dv:f:0.0694444  id:f:0.930556   cg:Z:4=1X2=1I38=1D5=1I5=1X13=

[user@cn3144]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. graphaligner.sh), which uses the input file 'graphaligner.in'. For example:



```

#!/bin/bash
module load graphaligner/1.0.16
GraphAligner -g $GRAPHALIGNER_TEST_DATA/graph.gfa -t $SLURM_CPUS_PER_TASK \
    -f $GRAPHALIGNER_TEST_DATA/read.fa -a aln.gaf -x vg

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=2 graphaligner.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. graphaligner.swarm). For example:



```

GraphAligner -g input1/graph.gfa -t $SLURM_CPUS_PER_TASK -f input1/read.fa -a output1/aln.gaf -x vg
GraphAligner -g input2/graph.gfa -t $SLURM_CPUS_PER_TASK -f input2/read.fa -a output2/aln.gaf -x vg

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f graphaligner.swarm -g 4 -t 2 --module graphaligner/1.0.16
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g #  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t #  Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module graphaligner  Loads the graphaligner module for each subjob in the swarm 
 | |
 | |
 | |








