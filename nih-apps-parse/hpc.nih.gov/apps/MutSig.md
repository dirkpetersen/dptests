

document.querySelector('title').textContent = 'MutSig on Biowulf';
MutSig on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
 |



MutSig analyzes lists of mutations discovered in DNA sequencing, to identify genes that were mutated more often than expected by chance given background mutation processes.




MutSigCV starts from the observation that the data is very sparse, and that there are usually too few silent mutations in a gene for its background mutation rate (BMR) to be estimated with any confidence. MutSigCV improves the BMR estimatation by pooling data from 'neighbor' genes in covariate space. These neighbor genes are chosen on the basis of having similar genomic properties to the central gene in question: properties such as DNA replication time, chromatin state (open/closed), and general level of transcription activity (e.g. highly transcribed vs. not transcribed at all). These genomic parameters have been observed to strongly correlate (co-vary) with background mutation rate. For instance, genes that replicate early in S-phase tend to have much lower mutation rates than late-replicating genes. Genes that are highly transcribed also tend to have lower mutation rates than unexpressed genes, due in part to the effects of transcription-coupled repair (TCR). Genes in closed chromatin (as measured by HiC or ChipSeq) have higher mutation rates than genes in open chromatin. Incorporating these covariates into the background model substantially reduces the number of false-positive findings.



### References:


* Lawrence, M. et al.
 [**Mutational heterogeneity in cancer and the search for new cancer-associated genes.**](https://www.ncbi.nlm.nih.gov/pubmed/23770567)
*Nature 499, 214-218 (2013).*


Documentation
* [MutSig Home](http://www.broadinstitute.org/cancer/cga/mutsig)


Important Notes
* Module Name: MutSig (see [the modules page](/apps/modules.html) for more information)
* Singlethreaded
* Environment variables set 
	+ MUTSIG\_HOME = path to the chosen version of MutSig
	+ MUTSIG\_EX = directory contain useful example files
	+ MUTSIG\_REF = directory contain requred reference files* Example files in $MUTSIG\_EX* Reference data in $MUTSIG\_REF



You will need to create a mutations file in MAF format.


Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive --mem=10g**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$
[user@cn3144 ~]$ **ml MutSig**
[user@cn3144 ~]$ **MutSigCV \
$MUTSIG\_EX/LUSC.mutations.maf \
$MUTSIG\_REF/exome\_full192.coverage.txt \
$MUTSIG\_REF/gene.covariates.txt \
output \
$MUTSIG\_REF/mutation\_type\_dictionary\_file.txt \
$MUTSIG\_REF/chr\_files\_hg19**
[user@cn3144 ~]$ **ls -l**
-rw-r--r-- 1 user user      510 Aug 29 12:47 output.categs.txt
-rw-r--r-- 1 user user  8989673 Aug 29 12:47 output.coverage.txt
-rw-r--r-- 1 user user 38760466 Aug 29 12:47 output.mutations.txt
-rw-r--r-- 1 user user      750 Aug 29 12:43 output.mutcateg_discovery.txt
-rw-r--r-- 1 user user  1397350 Aug 29 13:09 output.sig_genes.txt
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```

Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. MutSig.sh). For example:



```

#!/bin/bash
module load MutSig
MutSigCV \
  $MUTSIG_EX/LUSC.mutations.maf \
  $MUTSIG_REF/exome_full192.coverage.txt \
  $MUTSIG_REF/gene.covariates.txt \
  output \
  $MUTSIG_REF/mutation_type_dictionary_file.txt \
  $MUTSIG_REF/chr_files_hg19

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --mem=10g MutSig.sh
```


10 GB memory is sufficient for this example job. You may need to increase the memory allocation for your own MutSig jobs. 





