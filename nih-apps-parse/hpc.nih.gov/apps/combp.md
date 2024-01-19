

document.querySelector('title').textContent = 'combp on Biowulf';
combp on Biowulf


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



Combp is a library to combine, analyze, group and correct p-values in BED files. Unique tools involve correction for spatial autocorrelation. This is useful for ChIP-Seq probes and Tiling arrays, or any data with spatial correlation. 








### References:


* Pedersen BS, Schwartz DA, Yang IV, Kechris KJ. *Comb-p: software for combining, analyzing, grouping and correcting spatially correlated P-values.*
Bioinformatics. 2012 Nov 15;28(22):2986-8. doi: 10.1093/bioinformatics/bts545. Epub 2012 Sep 5.
 [PubMed](https://www.ncbi.nlm.nih.gov/pubmed/22954632) | 
 [Journal](https://academic.oup.com/bioinformatics/article/28/22/2986/240603)


Documentation
* combp Main Site:[Main Site](https://github.com/brentp/combined-pvalues)


Important Notes
* Module Name: combp (see [the modules page](/apps/modules.html) for more information)
 
```

	comb-p --help
```
* Environment variables set 
	+ $COMBP\_TESTDATA



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive --cpus-per-task=2 --mem=2G**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load combp**
[user@cn3144 ~]$ **comb-p**
Tools for viewing and adjusting p-values in BED files.

   Contact: Brent Pedersen - bpederse@gmail.com
   License: BSD

To run, indicate one of:

   acf       - calculate autocorrelation within BED file
   slk       - Stouffer-Liptak-Kechris correction of correlated p-values
   fdr       - Benjamini-Hochberg correction of p-values
   peaks     - find peaks in a BED file.
   region_p  - generate SLK p-values for a region (of p-values)
   filter    - filter region_p output on size and p and add coef/t
   hist      - plot a histogram of a column and check for uniformity.
   manhattan - a manhattan plot of values in a BED file.
   pipeline  - run acf, slk, fdr, peaks, region_p in succesion

NOTE: most of these assume *sorted* BED files.
SEE: https://github.com/brentp/combined-pvalues for documentation

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. combp.sh). For example:



```

#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G
#SBATCH --time=2:00:00
#SBATCH --partition=norm

set -e
module load combp
cp $COMBP_TESTDATA/pvals.bed .
comb-p acf -d 1:500:50 -c 5 pvals.bed > acf.txt

```

 Submit the job:

```
sbatch combp.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. job.swarm). For example:



```


       cd dir1; comb-p peaks --seed 0.05 --dist 1000 pvals1.bed >pvals.peaks1.bed
       cd dir2; comb-p peaks --seed 0.05 --dist 1000 pvals2.bed >pvals.peaks2.bed

    
```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f job.swarm [-g #] --module combp
```

where
 

|  |  |
| --- | --- |
| -g *#*  | Number of Gigabytes of memory required for each process
 (1 line in the swarm command file)  |
| --module  | Loads the module for each subjob in the swarm  |










