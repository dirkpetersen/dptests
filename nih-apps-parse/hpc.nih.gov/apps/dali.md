

document.querySelector('title').textContent = 'Dali on Biowulf';
Dali on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
[AlphaFold searching](#alphafold) 
 |



The three-dimensional co-ordinates of each protein are used to calculate residue - residue distance matrices. 



### References:


* Holm J.
 [**Using Dali for protein structure comparison.**](https://pubmed.ncbi.nlm.nih.gov/32006276/)
*Methods Mol. Biol. 2112, 29-42*


Documentation
* [Dali Main Site](http://ekhidna2.biocenter.helsinki.fi/dali/)


Important Notes
* Module Name: dali (see [the modules page](/apps/modules.html) for more information)
 * singlethreaded and MPI
* Environment variables set 
	+ DALI\_HOME
	+ DALI\_AF* Example files in $DALI\_HOME/example/* Reference data in /pdb/ and $DALI\_AF* Test script:  


```

#!/bin/bash
#SBATCH -J dali_test --ntasks=4 --nodes=1
rm -rf test
ml dali
$DALI_HOME/test.csh

```



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive --ntasks=4 --nodes=1**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load dali**
[user@cn3144 ~]$ **cp /pdb/pdb/pp/pdb1ppt.ent.gz .**
[user@cn3144 ~]$ **cp /pdb/pdb/bb/pdb1bba.ent.gz .**
[user@cn3144 ~]$ **import.pl --pdbfile pdb1ppt.ent.gz --pdbid 1ppt --dat ./**
[user@cn3144 ~]$ **import.pl --pdbfile pdb1bba.ent.gz --pdbid 1bba --dat ./**
[user@cn3144 ~]$ **dali.pl --pdbfile1 pdb1ppt.ent.gz --pdbfile2 pdb1bba.ent.gz --dat1 ./ --dat2 ./ --outfmt "summary,alignments"**
[user@cn3144 ~]$ **cat mol1A.txt**
# Job: test
# Query: mol1A
# No:  Chain   Z    rmsd lali nres  %id PDB  Description
   1:  mol2-A  3.6  1.8   33    36   39   MOLECULE: BOVINE PANCREATIC POLYPEPTIDE;

# Pairwise alignments

No 1: Query=mol1A Sbjct=mol2A Z-score=3.6

DSSP  LLLLLLLLLLLLLHHHHHHHHHHHHHHHHHHLLlll
Query GPSQPTYPGDDAPVEDLIRFYDNLQQYLNVVTRhry   36
ident  |  | |||| |  |        |  | |  ||
Sbjct APLEPEYPGDNATPEQMAQYAAELRRYINMLTRpry   36
DSSP  LLLLLLLLLLLLLLLHHHHHHHHHHHHHHHHLLlll

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. dali.sh). For example:



```

#!/bin/bash
module load dali
import.pl --pdbfile pdb1ppt.ent.gz --pdbid 1ppt --dat ./
import.pl --pdbfile pdb1bba.ent.gz --pdbid 1bba --dat ./
dali.pl --pdbfile1 pdb1ppt.ent.gz --pdbfile2 pdb1bba.ent.gz --dat1 ./ --dat2 ./ --outfmt "summary,alignments"

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch dali.sh
```

MPI batch job
In certain circumstances, dali can be accelerated using MPI. To do so, include **--np $SLURM\_NTASKS** with the command, and submit the job using **--ntasks=*#* --nodes=1** , where ***#*** is the number of MPI tasks requested. MPI only works on a single node, so # must be less than the maximum number of cpus available on a single node. At present the maximum is 128; however, most nodes have only 56 cpus and so jobs requesting more than 56 cpus may wait a considerable time in the queue.



```

...
dali.pl **--np $SLURM\_NTASKS** ...
...

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --ntasks=32 --nodes=1 dali.sh
```

AlphaFold searching
[Running with the AlphaFold database](http://ekhidna2.biocenter.helsinki.fi/dali/digest.html):



```

#!/bin/bash
module load dali

zcat /pdb/pdb/fd/pdb1fd3.ent.gz > 1fd3.pdb
import.pl --pdbfile 1fd3.pdb --pdbid 1fd3 --dat ./ --clean

dali.pl \
  --title "my search" \
  --cd1 1fd3B \
  --dat1 ./ \
  --db ${DALI_AF}/Digest/HUMAN.list \
  --BLAST_DB ${DALI_AF}/Digest/AF.fasta \
  --repset ${DALI_AF}/Digest/HUMAN_70.list \
  --dat2 ${DALI_AF}/DAT/ \
  --clean \
  --hierarchical \
  --oneway \
  --np ${SLURM_NTASKS}

```

Type **ls ${DALI\_AF}/Digest** to see all the lists.


**NOTES:**


* The import.pl process may fail due to formatting errors in the original pdb file. Make sure that it completes normally and individual .dat files are created for each chain. For the example (1fd3.pdb), there are four chains in the original pdb:
* The value of --cd1 must match the desired chain. So for the B chain of 1fd3.pdb, use 1fd3B
* Make sure to use ${SLURM\_NTASKS} for the value of --np; also be sure to allocate the job with --ntasks.








