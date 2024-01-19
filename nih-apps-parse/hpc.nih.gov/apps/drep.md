

document.querySelector('title').textContent = 'drep on Biowulf';
drep on Biowulf


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



dRep is a python program for rapidly comparing large numbers of genomes. dRep can also "de-replicate" a genome set by identifying groups of highly similar genomes and choosing the best representative genome for each genome set.






### References:


* Matthew R Olm et al.**dRep: a tool for fast and accurate genomic comparisons that enables improved genome recovery from metagenomes through de-replication**
 *The ISME Journal volume 11, pages2864–2868 (2017)* [Journal](https://www.nature.com/articles/ismej2017126)


Documentation
* drep Github:[Github](https://github.com/MrOlm/drep)


Important Notes
* drep can be run as:
 
```

	dRep --help
	
```
* Figures are located within the work directory under figures folder:
 
```

	ls out/figures/
	Clustering_scatterplots.pdf
	Cluster_scoring.pdf
	Primary_clustering_dendrogram.pdf
	Secondary_clustering_dendrograms.pdf
	Winning_genomes.pdf
	
```


![](/images/Primary_clustering_dendrogramK.png)




Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive --cpus-per-task=10 --mem=20G**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144]$ **module load drep**
[user@cn3144]$ **mkdir /data/$USER/drep\_test/**
[user@cn3144]$ **cd /data/$USER/drep\_test/**
[user@cn3144]$ **cp -r ${DREP\_TEST\_DATA:-none}/\* .**
[user@cn3144]$ **dRep --help**
...::: dRep v3.2.2 :::...

  Matt Olm. MIT License. Banfield Lab, UC Berkeley. 2017 (last updated 2020)

  See https://drep.readthedocs.io/en/latest/index.html for documentation
  Choose one of the operations below for more detailed help.

  Example: dRep dereplicate -h

  Commands:
    compare            -> Compare and cluster a set of genomes
    dereplicate        -> De-replicate a set of genomes
    check_dependencies -> Check which dependencies are properly installed


[user@cn3144]$ **dRep dereplicate out3 -g ./test/\*fasta**
***************************************************
    ..:: dRep dereplicate Step 1. Filter ::..
***************************************************

Will filter the genome list
2 genomes were input to dRep
Calculating genome info of genomes
100.00% of genomes passed length filtering
Running prodigal
Running checkM
GenomeInfo has no values over 1 for contamination- these should be 0-100, not 0-1!
100.00% of genomes passed checkM filtering
***************************************************
    ..:: dRep dereplicate Step 2. Cluster ::..
***************************************************

Running primary clustering
Running pair-wise MASH clustering
2 primary clusters made
Running secondary clustering
Running 2 ANImf comparisons- should take ~ 0.2 min
Step 4. Return output
***************************************************
    ..:: dRep dereplicate Step 3. Choose ::..
***************************************************

Loading work directory
GenomeInfo has no values over 1 for contamination- these should be 0-100, not 0-1!
***************************************************
    ..:: dRep dereplicate Step 4. Evaluate ::..
***************************************************

will provide warnings about clusters
0 warnings generated: saved to /spin1/home/linux/apptest1/out3/log/warnings.txt
will produce Widb (winner information db)
Winner database saved to /spin1/home/linux/apptest1/out3data_tables/Widb.csv
***************************************************
    ..:: dRep dereplicate Step 5. Analyze ::..
***************************************************

making plots 1, 2, 3, 4, 5, 6
Plotting primary dendrogram
Plotting secondary dendrograms
Plotting MDS plot
Plotting scatterplots
GenomeInfo has no values over 1 for contamination- these should be 0-100, not 0-1!
Plotting bin scorring plot
GenomeInfo has no values over 1 for contamination- these should be 0-100, not 0-1!
Plotting winning genomes plot...

$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    ..:: dRep dereplicate finished ::..

Dereplicated genomes................. /spin1/home/linux/apptest1/out3/dereplicated_genomes/
Dereplicated genomes information..... /spin1/home/linux/apptest1/out3/data_tables/Widb.csv
Figures.............................. /spin1/home/linux/apptest1/out3/figures/
Warnings............................. /spin1/home/linux/apptest1/out3/log/warnings.txt

$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

[user@cn3144]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. drep.sh). For example:



```


#!/bin/bash

#SBATCH --job-name=drep_run
#SBATCH --time=2:00:00

#SBATCH --partition=norm
#SBATCH --nodes=1
#SBATCH --mem=20g
#SBATCH --cpus-per-task=4

cd /data/$USER/drep_test
module load drep
dRep dereplicate out2 -g ./test/*fasta 




```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch drep.sh
```









