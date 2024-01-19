

document.querySelector('title').textContent = 'DANPOS on Biowulf';
DANPOS on Biowulf


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



A toolkit for Dynamic Analysis of Nucleosome and Protein Occupancy by Sequencing.




Features

 DANPOS provides different tools to analyze chromatin feature in different manner:
 * **dpos** analyzes changes in the location, fuzziness, and occupancy at each nucleosome or protein binding position. The functions is designed to analyze nucleosomes, but can be also very useful for analyzing other proteins whose binding pattern is similar to that of nucleosome, such as the protein MeCP2.
 * **dpeak** analyzes change in width, height and total signal of each enriched peak. A peak may contains multiple positions. This function is designed to analyze transcription factor and is also useful for defining some histone modification peaks.
 * **dregion** analyzes change in width, summit height, and total signal in each enriched region between samples. A region may contains multiple peaks, this function is designed to analyze chromatin features such as the histone modification H3K4me3.
 * **dtriple** is a combination of the Dpos, Dpeak, and Dregion functions. When there is no knowledge about characteristics of the sequencing target, this function provide users a way to try all the three algorithms above.
 * **profile** is a tool in DANPOS for analyzing the distribution of a chromatin feature flanking each given group of genomic sites or regions.
 * **stat** is a tool in DANPOS to do statistical analysis for positions, peaks, or regions.
 * **wiq** is a tool in DANPOS to do genome wide quantile normalization for data in wiggle format file.
 * **wig2wiq** is a tool for converting .wig format file to .wiq format.




### References:


* Kaifu Chen, etc. *DANPOS: Dynamic Analysis of Nucleosome Position and Occupancy by Sequencing.*Genome Research. 2012. doi:10.1101/gr.142067.112 
 [PubMed](https://www.ncbi.nlm.nih.gov/pubmed/23193179) | 
 [PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3561875/) | 
 [Journal](https://genome.cshlp.org/content/23/2/341.long)


Documentation
* DANPOS Main Site:[Main Site](https://sites.google.com/site/danposdoc/news)
* DANPOS Github:[Github](https://github.com/sklasfeld/DANPOS3)


Important Notes
* Module Name: DANPOS (see [the modules page](/apps/modules.html) for more information)
 * Current DANPOS command lines could be run as:
 
```

	danpos.py -h

```



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

[user@cn3144 ~]$ **module load DANPOS**
[user@cn3144 ~]$ **danpos.py --help**
danpos version 3.0.0
For help information for each function, try:
python danpos.py  -h

Functions:
 dpos:
 analyze each protein-binding position (~100
 to ~200bp wide) across the whole genome,
 e.g. nucleosome positions.
 dpeak:
 analyze each protein-binding peak (~1 to ~1kb
 wide) across the whole genome, e.g. protein
 that binds accruately to some specific motifs.
 dregion:
 analyze each protein-binding region (~1 to
 ~10kb wide) across the whole genome, e.g.
 some histone modifications.
 dtriple:
 Do analysis at all three levels including each
 region, peak, and position. Would be useful
 when little is known about the potential binding
 pattern.
 profile:
 analyze wiggle format occupancy or differential
 signal profile relative to gene structures or
 bed format genomic regions.
 wiq:
 normalize wiggle files to have the same quantiles (Quantile normalization).
 wig2wiq:
 convert wiggle format file to wiq format.
 stat:
 some statistics for positions, peaks, or regions.
 selector:
 select a subset of positions, peaks, or regions
 by value ranges or gene structures neighboring.
 valuesAtRanks:
 retrieve position, peak, or region values by ranks.

Kaifu Chen, et al. chenkaifu@gmail.com, Li lab, Biostatistics department, Dan L. Duncan cancer center, Baylor College of Medicine.
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. DANPOS.sh). For example:




hljs.highlightAll();

```


#!/bin/bash
set -e
module load DANPOS
danpos.py dops sampleA


```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=2 --mem=2g DANPOS.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. DANPOS.swarm). For example:



```

cd dir1;danpos.py dpos sampleA
cd dir2;danpos.py dpos sampleB

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f DANPOS.swarm [-g #] [-t #] --module DANPOS
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module DANPOS Loads the DANPOS module for each subjob in the swarm
 | |
 | |
 | |








