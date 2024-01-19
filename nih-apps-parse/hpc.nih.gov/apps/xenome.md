

document.querySelector('title').textContent = 'xenome: a tool for classifying reads from xenograft sources';
**xenome: a tool for classifying reads from xenograft sources**


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



Xenome performs fast, accurate and specific classification
of xenograft-derived sequence read data. It handles a mixture of 
reads arising from the host and reads
arising from the graft and separates the
two, thus allowing for more precise analysis to be performed.



### References:


* T.Conway, J.Wazny, A.Bromage, M.Tymms, D.Sooraj, E.D.Williams and B.Beresford-Smith   

*Xenome—a tool for classifying reads from xenograft samples.*   

[Bioinformatics](https://academic.oup.com/bioinformatics/article/28/12/i172/269972)  2012, **28**, i172–i178. DOI: 10.1093/bioinformatics/bts236


Documentation
* [xenome user manual](https://github.com/data61/gossamer/blob/master/docs/xenome.md)


Important Notes
* Module Name: xenome (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set
	+ **XENOME\_HOME**  xenome installation directory
	+ **XENOME\_BIN**       xenome executable directory
	+ **XENOME\_DATA**    xenome sample data directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=45g --cpus-per-task=8**
[user@cn3316 ~]$**module load xenome**

```

Xenome has two distinct stages, which are embodied in two separate commands: 'index' and 'classify'. Before reads can be classified, an index must be constructed from the graft and host reference sequences. The references must be in FASTA format, and may optionally be compressed (gzip).   
  

1) Build a xenome index from the host (H, mouse) and graft (G, human) genome reference sequences:

```

[user@cn3316 ~]$**ln -s /fdb/KNIFE/hg19/hg19\_genome.fa human.fa**
[user@cn3316 ~]$**ln -s /fdb/muTect/mm9.fa mouse.fa**
[user@cn3316 ~]$**xenome index -M 24 -T 8 -P idx -H mouse.fa -G human.fa -K 20 -M 45**
581200	1416552	4065928400	0.0348396
970281	2431208	4065928400	0.0597947
1121297	2839948	4065928400	0.0698475
1395390	3545570	4065928400	0.087202
2151142	5956366	4065928400	0.146495
...
1013843220	4063007338	4065928400	99.9282
1014221284	4064055913	4065928400	99.9539
1014677858	4065104499	4065928400	99.9797
1014924662	4065928400	4065928400	100

```

This will create a number of related files which can be identified by a user-specified prefix, e.g. 'idx' in the above command.   
   

**NOTE: step 1) will take several hours to complete. For the particular case of mouse as a host and human a as graft, this step was pre-computed already, so instead of re-building the index one can simply download the results into the current folder:**

```

[user@cn3316 ~]$**cp $XENOME\_DATA/\* .**

```
 
2) Once an index is available, reads can be classified according to whether they appear to contain graft or host material. In the simplest case, xenome can classify each read from a single source file individually.

```

[user@cn3316 ~]$**ln -s $XENOME\_DATA/SRR4254643\_mouse.fastq** 
[user@cn3316 ~]$**xenome classify -P idx -i SRR4254643\_mouse.fastq** 
Statistics
B	G	H	M	count	percent	class
0	0	0	0	70719	0.184239	"neither"
0	0	0	1	14071	0.0366581	"both"
0	0	1	0	22611	0.0589066	"definitely host"
0	0	1	1	438587	1.14262	"probably host"
0	1	0	0	58415	0.152184	"definitely graft"
0	1	0	1	8877415	23.1276	"probably graft"
0	1	1	0	6469	0.0168532	"ambiguous"
0	1	1	1	174335	0.454181	"ambiguous"
1	0	0	0	77780	0.202634	"both"
1	0	0	1	19153	0.0498978	"probably both"
1	0	1	0	6430	0.0167516	"definitely host"
1	0	1	1	2324237	6.05515	"probably host"
1	1	0	0	75427	0.196504	"definitely graft"
1	1	0	1	25756418	67.1012	"probably graft"
1	1	1	0	1855	0.00483268	"ambiguous"
1	1	1	1	460543	1.19982	"ambiguous"

Summary
count	percent	class
34767675	90.5775	graft
2791865	7.27342	host
111004	0.28919	both
70719	0.184239	neither
643202	1.67568	ambiguous


```

The latter command will also create files: 

```

ambiguous.fastq  graft.fastq  neither.fastq
both.fastq       host.fastq

```

End the interactive session:

```

[user@cnR3316 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```

Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. xenome.sh). For example:



```

#!/bin/bash
module load xenome
cp $XENOME_DATA/* .
ln -s /fdb/muTect/mm9.fa             mouse.fa
ln -s /fdb/KNIFE/hg19/hg19_genome.fa human.fa
xenome xenome -P idx -i SRR4254643_mouse.fastq

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] xenome.sh
```





