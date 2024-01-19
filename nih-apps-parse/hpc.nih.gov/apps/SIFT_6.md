

document.querySelector('title').textContent = 'SIFT (6.x) on Biowulf';
SIFT (6.x) on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
 |



SIFT predicts whether an amino acid substitution affects protein function based on sequence homology and the physical properties of amino acids. SIFT can be applied to naturally occurring nonsynonymous polymorphisms and laboratory-induced missense mutations.



### References:


* Sim NL, Kumar P, Hu J, Henikoff S, Schneider G, Ng PC.
	[SIFT web server: predicting effects of amino acid substitutions on proteins.](https://www.ncbi.nlm.nih.gov/pubmed/22689647)*Nucleic Acids Res. 2012 Jul;40*


Documentation
The documentation here pertains to SIFT version 6.x. For version 5.x, [click here](SIFT_5.html).


* [SIFT website](http://sift.bii.a-star.edu.sg/)


Important Notes
* Module Name: SIFT (see [the modules page](/apps/modules.html) for more information)
* Multithreaded
* Environment variables set 
	+ **BLIMPS\_DIR**: directory where blimps is installed
	 + **NCBI**: directory where NCBI blast is installed
	 + **SIFTDB**: directory where SIFT reference files are held
	 + **SIFTHOME**: directory where SIFT is installed
	 + **SIFT\_SCRATCHDIR**: directory where SIFT output is written* Example files in $SIFT/home* Reference data in /fdb/SIFT/


SIFT can require a large amount of disk space. The environment variable **$SIFT\_SCRATCHDIR** is set to
**/lscratch/$SLURM\_JOB\_ID** by default, but can be changed.



```
export SIFT_SCRATCHDIR=/path/to/new/tmp/area
```

SIFT can use the following databases for protein alignment:


* /fdb/blastdb/nr -- NCBI non-redundant
* $SIFTDB/uniref90.fa -- UniRef 90
* $SIFTDB/uniprotkb\_swissprot.fa -- SwissProt


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

[user@cn3144 ~]$ module load SIFT/6.2.1
[user@cn3144 ~]$ cp $SIFTHOME/test/lacI.fasta .
[user@cn3144 ~]$ SIFT_for_submitting_fasta_seq.csh lacI.fasta $SIFTDB/uniref90.fa -
[user@cn3144 ~]$ mv $SIFT_SCRATCHDIR/lacI.* .

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. SIFT.sh). For example:



```

#!/bin/bash
module load SIFT/5.2.2
cp $SIFTHOME/test/lacI.fasta .
SIFT_for_submitting_fasta_seq.csh lacI.fasta $SIFTDB/uniref90.fa -
mv $SIFT_SCRATCHDIR/lacI.* .

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] TEMPLATE.sh
```





