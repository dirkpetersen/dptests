

document.querySelector('title').textContent = 'tRNAscan-SE on Biowulf';
tRNAscan-SE on Biowulf


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



 tRNAscan-SE 2.0 has advanced the state-of-the-art methodology in tRNA gene detection and functional prediction, captured by rich new content of the companion Genomic tRNA Database


### References:


* Lowe, T.M. and Chan, P.P. (2016) tRNAscan-SE On-line: Search and Contextual Analysis of Transfer RNA Genes. Nucl. Acids Res. 44: W54-57.


Documentation
* [tRNAscan-SE Main Site](http://lowelab.ucsc.edu/tRNAscan-SE/)


Important Notes
* Module Name: trnascan-se (see [the modules page](/apps/modules.html) for more information)
* Singlethreaded app
* Example files in /usr/local/apps/trnascan-se/TEST\_DATA



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

[user@cn3144 ~]$ **module load trnascna-se**
[+] Loading trnascan-se  2.0.0  on cn0166 

[user@cn3144 ~]$ **cp /usr/local/apps/trnascan-se/TEST\_DATA/ExampleSequences.fa .**

[user@cn3144 ~]$ **tRNAscan-SE ExampleSequences.fa**

tRNAscan-SE v.2.0 (December 2017) - scan sequences for transfer RNAs
Copyright (C) 2017 Patricia Chan and Todd Lowe
                   University of California Santa Cruz
Freely distributed under the GNU General Public License (GPLv3)

------------------------------------------------------------
Sequence file(s) to search:        ExampleSequences.fa
Search Mode:                       Eukaryotic
Results written to:                Standard output
Output format:                     Tabular
Searching with:                    Infernal First Pass->Infernal
Isotype-specific model scan:       Yes
Covariance model:                  /usr/local/apps/trnascan-se/2.0.0/lib/tRNAscan-SE/models/TRNAinf-euk.cm
                                   /usr/local/apps/trnascan-se/2.0.0/lib/tRNAscan-SE/models/TRNAinf-euk-SeC.cm
Infernal first pass cutoff score:  10

Temporary directory:               /tmp
------------------------------------------------------------

Status: Phase I: Searching for tRNAs with HMM-enabled Infernal
Status: Phase II: Infernal verification of candidate tRNAs detected with first-pass scan
Sequence		tRNA	Bounds	tRNA	Anti	Intron Bounds	Inf	      
Name    	tRNA #	Begin	End	Type	Codon	Begin	End	Score	Note
--------	------	-----	------	----	-----	-----	----	------	------
MySeq1  	1	13 	85 	Thr	TGT	0	0	78.0	
MySeq2  	1	6  	79 	Arg	TCT	0	0	75.1	
MySeq3  	1	14 	114	Ser	CGA	51	69	71.8	
MySeq4  	1	6  	88 	Leu	AAG	0	0	65.0	
MySeq5  	1	3  	89 	SeC	TCA	0	0	146.9	
MySeq6  	1	7  	92 	Lys	CTT	0	0	72.1	

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. trnascan-se.sh). For example:



```

#!/bin/bash
set -e
module load trnascan-se
tRNAscan-SE ExampleSequences.fa

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch trnascan-se.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. trnascan-se.swarm). For example:



```

tRNAscan-SE ExampleSequences1.fa > out1.txt
tRNAscan-SE ExampleSequences2.fa > out2.txt
tRNAscan-SE ExampleSequences3.fa > out3.txt

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f trnascan-se.swarm --module trnascan-se
```

where


|  |  |
| --- | --- |
| --module trnascan-se Loads the trnascan-se module for each subjob in the swarm 
 | |








