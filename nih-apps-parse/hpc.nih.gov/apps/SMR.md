

document.querySelector('title').textContent = 'SMR: Summary-data-based Mendelian Randomization tool on Biowulf';
**SMR: Summary-data-based Mendelian Randomization tool on Biowulf**


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



SMR (Summary-based Mendelian Randomization) software 
integrates summary-level data from genome-wide association studies (GWAS) 
with data from expression quantitative trait locus (eQTL) studies 
to identify genes whose expression levels are associated with a complex trait because of pleiotropy. 
It implements the SMR & HEIDI methods to test for pleiotropic association 
between the expression level of a gene and a complex trait of interest 
using summary-level data from GWAS and expression quantitative trait loci (eQTL) studies (Zhu et al. 2016 Nat Genet). 
The methodology can be interpreted as an analysis to test if the effect size of a SNP 
on the phenotype is mediated by gene expression. 
This tool can therefore be used to prioritize genes underlying GWAS hits for follow-up functional studies.



### References:


* Zhihong Zhu, Futao Zhang, Han Hu, Andrew Bakshi, Matthew R Robinson, Joseph E Powell,
Grant W Montgomery, Michael E Goddard, Naomi R Wray, Peter M Visscher and Jian Yang,
"Integration of summary data from GWAS and eQTL studies predicts complex trait gene targets", Nature Genetics 2016, v.48(5), 481-489


Documentation
* SMR Main Site: [SMR](http://cnsgenomics.com/software/smr)


Important Notes
* Module Name: SMR (see [the modules page](/apps/modules.html) for more information)
* Multithreaded 
* Unusual environment variables set
	+ **SMR\_DIR** SMR installation directory* Example files in **$SMR\_DIR**/TEST\_DATA



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

[user@cn3144 ~]$ **module load SMR**
[+] Loading SMR  1.03
[user@cn3144 ~]$ **smr --beqtl-summary $SMR\_DIR/TEST\_DATA/westra\_eqtl\_hg18 --query 5.0e-8 --snp rs123**
*******************************************************************
* Summary-data-based Mendelian Randomization (SMR)
* version 1.03
* (C) 2015 Futao Zhang, Zhihong Zhu and Jian Yang
* The University of Queensland
* MIT License
*******************************************************************
Analysis started: 13:27:52,Tue Sep 29,2020

Options:
--beqtl-summary /usr/local/apps/SMR/TEST_DATA/westra_eqtl_hg18
--query   5.00e-08
--snp rs123


Reading eQTL summary data...
Reading eQTL probe information from [/usr/local/apps/SMR/TEST_DATA/westra_eqtl_hg18.epi].
5967 Probes to be included from [/usr/local/apps/SMR/TEST_DATA/westra_eqtl_hg18.epi].
Reading eQTL SNP information from [/usr/local/apps/SMR/TEST_DATA/westra_eqtl_hg18.esi].
506193 SNPs to be included from [/usr/local/apps/SMR/TEST_DATA/westra_eqtl_hg18.esi].
rs123 is extracted.
Reading eQTL summary data from [/usr/local/apps/SMR/TEST_DATA/westra_eqtl_hg18.besd].
This is an old file format. Please use --make-besd to update the file format.
eQTL summary data of 5967 Probes and 1 SNPs to be included from [/usr/local/apps/SMR/TEST_DATA/westra_eqtl_hg18.besd].
Extracted results of 1 SNPs have been saved in the file [smr.txt].

Analysis completed: 13:27:53,Tue Sep 29,2020
Computational time: 0:0:1
[user@cn3144 ~]$ **cat smr.txt**
SNP     Chr     BP      A1      A2      Freq    Probe   Probe_Chr       Probe_bp        Gene    Orientation b       SE      p
rs123   7       24932971        A       C       NA      ILMN_1670145    7       24704613        DFNA5       -       0.117336        0.0193413       1.30632e-09
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. SMR.sh). For example:



```

#!/bin/bash
module load SMR       
smr --beqtl-summary $SMR_DIR/TEST_DATA/westra_eqtl_hg18 --query 5.0e-8 --snp rs123

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] SMR.sh
```







