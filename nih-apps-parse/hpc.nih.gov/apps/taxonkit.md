

document.querySelector('title').textContent = 'TaxonKit on Biowulf';
TaxonKit on Biowulf


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



A Cross-platform and Efficient NCBI Taxonomy Toolkit




Subcommands
* list: list taxon tree of given taxids
 * lineage: query lineage of given taxids (supporting STDIN)
 * reformat: reformat lineage (supporting STDIN)
 * name2taxid: query taxid by taxon scientific name (supporting STDIN)
 * taxid-changelog: create taxid changelog from dump archives




### References:


* Shen W, Xiong J. *TaxonKit: a cross-platform and efficient NCBI taxonomy toolkit*. doi:https://doi.org/10.1101/513523


Documentation
* taxonkit Main Site:[Main Site](https://bioinf.shenwei.me/taxonkit/)


Important Notes
* Module Name: taxonkit (see [the modules page](/apps/modules.html) for more information)



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

[user@cn3144 ~]$ **module load taxonkit**
[user@cn3144 ~]$ **taxonkit list --show-rank --show-name --ids 9605**
9605 [genus] Homo
  9606 [species] Homo sapiens
    63221 [subspecies] Homo sapiens neanderthalensis
    741158 [subspecies] Homo sapiens subsp. 'Denisova'
    2665952 [no rank] environmental samples
      2665953 [species] Homo sapiens environmental sample
  1425170 [species] Homo heidelbergensis

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. taxonkit.sh). For example:



```

#!/bin/bash
set -e
module load taxonkit
#This will download the whole taxon tree
taxonkit list --show-rank --show-name --ids 1 >t1.txt

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch  --mem=2g taxonkit.sh
```







