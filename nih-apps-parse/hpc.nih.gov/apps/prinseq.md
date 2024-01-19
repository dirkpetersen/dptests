

document.querySelector('title').textContent = 'PRINSEQ';
PRINSEQ


|  |
| --- |
| 
Quick Links
[On Helix](#helix)
[Interactive job on Biowulf](#int)
[Batch job on Biowulf](#sbatch)
[Swarm of jobs on Biowulf](#swarm)
[Documentation](#doc)
 |



PRINSEQ is a tool that generates summary statistics of sequence and quality data and that is used to filter, reformat and trim next-generation sequence data. It is particular designed for 454/Roche data, but can also be used for other types of sequence data.




On our systems, only the main program prinseq-lite.pl is set up for use. Please contact the HPC staff if you require more of PRINSEQ's functionality.



### References:


* Schmieder R and Edwards R: Quality control and preprocessing of metagenomic datasets. Bioinformatics 2011, 27:863-864.


There may be multiple versions of prinseq available. An easy way of selecting the version is to use [modules](/apps/modules.html). To see the modules available, type



```
module avail prinseq
```

To select a module, type



```
module load prinseq/[ver]
```

where [ver] is the version of choice.


### Environment variables set:


* $PATH
* $PRINSEQ\_HOME


On Helix

Sample session:

```
module load prinseq
prinseq-lite.pl -verbose -fastq $PRINSEQ_HOME/example/example1.fastq -ns_max_n 0 -out_good test_no_ns -out_bad test_with_ns

```

Interactive job on Biowulf

See the [Biowulf user guide for interactive jobs](https://hpc.nih.gov/docs/userguide.html#int).



Batch job on Biowulf
Create a batch input file (e.g. prinseq.sh), which uses the input file 'prinseq.in'. For example:



```
#!/bin/bash
module load prinseq
prinseq-lite.pl -verbose -fastq $PRINSEQ_HOME/example/example1.fastq -ns_max_n 0 -out_good test_no_ns -out_bad test_with_ns

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=1 prinseq.sh
```

Swarm of Jobs on Biowulf

Create a swarmfile following the [swarm guide](https://hpc.nih.gov/apps/swarm.html) using the example commands on this page.



Documentation
* PRINSEQ Main Site: <http://prinseq.sourceforge.net>






