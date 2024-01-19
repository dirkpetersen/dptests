

document.querySelector('title').textContent = 'snpEff on Biowulf';
snpEff on Biowulf


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


snpEff is a variant annotation and effect prediction tool. It annotates and predicts the effects of variants on genes (such as amino acid changes).


Typical usage :


* **Input:** The inputs are predicted variants (SNPs, insertions, deletions and MNPs). The input file is usually obtained as a result of a sequencing experiment, and it is usually in variant call format (VCF).
* **Output:** SnpEff analyzes the input variants. It annotates the variants and calculates the effects they produce on known genes (e.g. amino acid changes). A list of effects and annotations that SnpEff can calculate can be found here.


### References:


* [A program for annotating and predicting the effects of single nucleotide polymorphisms, SnpEff: SNPs in the genome of Drosophila melanogaster strain w1118; iso-2; iso-3.", Cingolani P, Platts A, Wang le L, Coon M, Nguyen T, Wang L, Land SJ, Lu X, Ruden DM. Fly (Austin). 2012 Apr-Jun;6(2):80-92. PMID: 22728672](https://www.ncbi.nlm.nih.gov/pubmed/?term=22728672)


Documentation
* [snpEff Homepage](http://snpeff.sourceforge.net/)* [snpSift Homepage](http://snpeff.sourceforge.net/SnpSift.html)* [snpEff Usage examples](http://snpeff.sourceforge.net/protocol.html)


Important Notes
* Module Name: snpEff (see [the modules page](/apps/modules.html) for more information)
* Single-threaded, java-dependent
* environment variables set 
	+ **$SNPEFF\_HOME** -- path to the installation directory
	+ **$SNPEFF\_JAR** -- the snpEff.jar file
	+ **$SNPEFF\_JARPATH** -- same as SNPEFF\_HOME
	+ **$SNPSIFT\_JAR** -- the SnpSift.jar file* Example files in $SNPEFF\_HOME/../protocols* Reference data in /fdb/dbNSFP2/


snpEff is a java application. See [https://hpc.nih.gov/development/java.html](/development/java.html) for information about running java applications.


To see the help menu, type



```
java -jar $SNPEFF_JAR
```

at the prompt.


By default, snpEff uses 1gb of memory. For large VCF input files, this may not be enough.
 To allocate 20gb of memory, use:



```
java -Xmx20g -jar $SNPEFF_JAR
```

In more detail,



```
java -Xmx20g -jar $SNPEFF_JAR -v [ database ] [ vcf file ] 
```

To see what databases are available, type:



```
ls $SNPEFF_HOME/data
```

### snpSift


SnpSift is a collection of tools to manipulate VCF (variant call format) files. Here's what you can do:


* **Filter:** You can filter using arbitrary expressions, for instance "(QUAL > 30) | (exists INDEL) | ( countHet() > 2 )". The actual expressions can be quite complex, so it allows for a lot of flexibility.
* **Annotate:** You can add 'ID' from another database (e.g. variants from dbSnp)
* **CaseControl:** You can compare how many variants are in 'case' and in 'control' groups. Also calculates p-values (Fisher exact test).
* **Intervals:** Filter variants that intersect with intervals.
* **Intervals (intidx):** Filter variants that intersect with intervals. Index the VCF file using memory mapped I/O to speed up the search. This is intended for huge VCF files and a small number of intervals to retrieve.
* **Join:** Join by generic genomic regions (intersecting or closest).
* **RmRefGen:** Remove reference genotype (i.e. replace '0/0' genotypes by '.')
* **TsTv:** Calculate transiton to transversion ratio.
* **Extract fields:** Extract fields from a VCF file to a TXT (tab separated) format.
* **Variant type:** Adds SNP/MNP/INS/DEL to info field. It also adds "HOM/HET" if there is only one sample.
* **GWAS Catalog:** Annotate using GWAS Catalog.
* **dbNSFP:** Annotate using dbNSFP: The dbNSFP is an integrated database of functional predictions from multiple algorithms (SIFT, Polyphen2, LRT and MutationTaster, PhyloP and GERP++, etc.)


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

[user@cn3144 ~]$ module load snpEff
[user@cn3144 ~]$ ln -s $SNPEFF_HOME/../protocols .
[user@cn3144 ~]$ java -Xmx12g -jar $SNPEFF_JAR -v -lof -motif -hgvs -nextProt GRCh37.71 protocols/ex1.vcf > ex1.eff.vcf
[user@cn3144 ~]$ cat ex.eff.vcf | java -jar $SNPSIFT_JAR filter "(Cases[0] = 3) & (Controls[0] = 0) & ((EFF[*].IMPACT = 'HIGH') | (EFF[*].IMPACT = 'MODERATE'))"  > ex1.filtered.vcf
[user@cn3144 ~]$ java -jar $SNPSIFT_JAR dbnsfp -v -db /fdb/dbNSFP2/dbNSFP2.9.txt.gz ex1.eff.vcf > file.annotated.vcf

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. snpEff.sh). For example:



```
#!/bin/bash
# -- this file is snpEff.sh --

module load snpEff
ln -s $SNPEFF_HOME/example/file.vcf .
java -Xmx${SLURM_MEM_PER_NODE}m -jar $SNPEFF_JAR -v hg19 file.vcf > file.eff.vcf
cat file.eff.vcf | java -jar $SNPSIFT_JAR filter "( EFF[*].IMPACT = 'HIGH' )" > file.filtered.vcf
java -jar $SNPSIFT_JAR dbnsfp -v -db /fdb/dbNSFP2/dbNSFP3.2a.txt.gz file.eff.vcf > file.annotated.vcf

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] snpEff.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. snpEff.swarm). For example:



```

java -Xmx${SLURM_MEM_PER_NODE}m -jar $SNPEFF_JAR -v hg19 file1.vcf > file1.eff.vcf
java -Xmx${SLURM_MEM_PER_NODE}m -jar $SNPEFF_JAR -v hg19 file2.vcf > file2.eff.vcf
java -Xmx${SLURM_MEM_PER_NODE}m -jar $SNPEFF_JAR -v hg19 file3.vcf > file3.eff.vcf
java -Xmx${SLURM_MEM_PER_NODE}m -jar $SNPEFF_JAR -v hg19 file4.vcf > file4.eff.vcf

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f snpEff.swarm [-g #] --module snpEff
```

where


|  |  |  |  |
| --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | --module snpEff Loads the snpEff module for each subjob in the swarm 
 | |
 | |






