

document.querySelector('title').textContent = 'sumtrees on Biowulf';
sumtrees on Biowulf


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



The sumtrees program summarize non-parameteric bootstrap or Bayesian posterior probability support for splits or clades on phylogenetic trees. 






### References:


* Sukumaran, J. and Mark T. Holder. *DendroPy: A Python library for phylogenetic computing.*Bioinformatics.2010.26: 1569-1571. [PubMed](https://www.ncbi.nlm.nih.gov/pubmed/20421198) | 
 [Journal](https://academic.oup.com/bioinformatics/article/26/12/1569/287181)


Documentation
* sumtrees Mainsite:[Mainsite](https://dendropy.org/programs/sumtrees.html)


Important Notes
* Module Name: sumtrees (see [the modules page](/apps/modules.html) for more information)
 * Current sumtrees command lines could be run as:
 
```

	sumtrees.py
	
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

[user@cn3144 ~]$ **module load sumtrees**
[user@cn3144 ~]$ **sumtrees.py --help**
[user@cn3144 ~]$ **sumtrees.py --min-clade-freq=0.95 --burnin=200 --output-tree-filepath=result.tre treefile1.tre treefile2.tre treefile3.tre**
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. sumtrees.sh). For example:




hljs.highlightAll();

```


#!/bin/bash
set -e
module load sumtrees
sumtrees.py -m4 -f0.95 -b200 -o result.tre treefile1.tre treefile2.tre treefile3.tre treefile4.tre



```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=4 --mem=2g sumtrees.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. sumtrees.swarm). For example:



```

cd dir1;sumtrees.py -f0.95 -b200 -o result.tre treefile1.tre treefile2.tre treefile3.tre
cd dir2;sumtrees.py --summary-target=mcct --burnin=200 --support-as-labels --output-tree-filepath=result2.tre treefile1.tre treefile2.tre treefile3.tre

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f sumtrees.swarm [-t #] [-g #] --module sumtrees
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module sumtrees Loads the sumtrees module for each subjob in the swarm
 | |
 | |
 | |










