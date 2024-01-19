

document.querySelector('title').textContent = 'Bedops on Biowulf';
Bedops on Biowulf


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



Bedops is a software package for manipulating and analyzing genomic interval
data. It contains tools to


* perform boolean and set operations on genome intervals
* carry out statistical calculations on genome intervals
* sort bed files consistently (important for a lot of operations)
* compress genomic data using it's own *starch* format
* conversions of various formats (bam, gff/gtf, wig, vcf) 
 to *bed* and *starch*


Operations can be parallelized efficiently by chromosome.


Each tool is designed to use unix input and output streams for building
efficient pipelines.





### References:


Shane Neph, et al. *BEDOPS: high-performance genomic feature operations*
 Bioinformatics 2012, 28:1919-1920.
 [Pubmed](http://www.ncbi.nlm.nih.gov/pubmed/22576172) | 
 [PMC](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3389768/) | 
 [Journal](http://bioinformatics.oxfordjournals.org/content/28/14/1919.long)

Documentation
* [GitHub](https://github.com/bedops/bedops)
* [Manual](http://bedops.readthedocs.org/en/latest/)
* [Overview](http://bedops.readthedocs.org/en/latest/content/overview.html#overview)
* [Examples](http://bedops.readthedocs.org/en/latest/content/usage-examples.html)



Important Notes
* Module Name: bedops (see [the modules 
 page](/apps/modules.html) for more information)





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

[user@cn3144 ~]$ **convert2bed --input=bam --output=starch --starch-note "wt/H3ac/a13" < bam\_inputfile > starch/a13.starch**
[user@cn3144 ~]$ **unstarch --note starch/a13.starch wt/H3ac/a13**
[user@cn3144 ~]$ **unstarch --elements starch/a13.starch**
[user@cn3144 ~]$ **unstarch --bases-uniq starch/a13.starch**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```




Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. bedops.sh). For example:



```

#!/bin/bash
set -e
module load bedops
sort-bed in.bed > out.bed
```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--mem=#] bedops.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. bedops.swarm). For example:



```

cd dir1; sort-bed in.bed > out.bed
cd dir2; sort-bed in.bed > out.bed
cd dir3; sort-bed in.bed > out.bed
cd dir4; sort-bed in.bed > out.bed


```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f bedops.swarm [-g #] --module bedops
```

where
 

|  |  |
| --- | --- |
| -g *#*  | Number of Gigabytes of memory required for each process 
 (1 line in the swarm command file)  |
| --module  | Loads the module for each subjob in the swarm  |




