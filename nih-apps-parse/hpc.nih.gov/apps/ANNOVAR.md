

document.querySelector('title').textContent = 'ANNOVAR on Biowulf';
ANNOVAR on Biowulf


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



ANNOVAR is an efficient software tool to utilize update-to-date information to functionally annotate genetic variants detected from diverse genomes (including human genome {hg18,hg19,hg38} as well as mouse, worm, fly, yeast and many others).



### References:


* [ANNOVAR: Functional annotation of genetic variants from next-generation sequencing data, Nucleic Acids Research, 38:e164, 2010)](http://www.ncbi.nlm.nih.gov/pubmed/?term=20601685).


Documentation
* [ANNOVAR Home](http://annovar.openbioinformatics.org/)
* [Wang K, Li M, Hakonarson H. ANNOVAR: Functional annotation of genetic variants from next-generation sequencing data Nucleic Acids Research, 38:e164, 2010](http://www.ncbi.nlm.nih.gov/pubmed/20601685)


Important Notes
* Module Name: annovar (see [the modules page](/apps/modules.html) for more information)
 * Multithreaded (table\_annovar.pl can utilize multiple threads)
 * Environment variables set 
	+ **ANNOVAR\_HOME**
	+ **ANNOVAR\_DATA*** Example files in **$ANNOVAR\_HOME/example*** Reference data in **/fdb/annovar/current/**


ANNOVAR takes text-based input files, where each line corresponds to one variant. On each line, the first five space- or tab- delimited columns represent chromosome, start position, end position, the reference nucleotides and the observed nucleotides.

Here is the example file **$ANNOVAR\_HOME/example/ex1.avinput**



```

1	948921	948921	T	C	comments: rs15842, a SNP in 5' UTR of ISG15
1	1404001	1404001	G	T	comments: rs149123833, a SNP in 3' UTR of ATAD3C
1	5935162	5935162	A	T	comments: rs1287637, a splice site variant in NPHP4
1	162736463	162736463	C	T	comments: rs1000050, a SNP in Illumina SNP arrays
1	84875173	84875173	C	T	comments: rs6576700 or SNP_A-1780419, a SNP in Affymetrix SNP arrays
1	13211293	13211294	TC	-	comments: rs59770105, a 2-bp deletion
1	11403596	11403596	-	AT	comments: rs35561142, a 2-bp insertion
1	105492231	105492231	A	ATAAA	comments: rs10552169, a block substitution
1	67705958	67705958	G	A	comments: rs11209026 (R381Q), a SNP in IL23R associated with Crohn's disease
2	234183368	234183368	A	G	comments: rs2241880 (T300A), a SNP in the ATG16L1 associated with Crohn's disease
16	50745926	50745926	C	T	comments: rs2066844 (R702W), a non-synonymous SNP in NOD2
16	50756540	50756540	G	C	comments: rs2066845 (G908R), a non-synonymous SNP in NOD2
16	50763778	50763778	-	C	comments: rs2066847 (c.3016_3017insC), a frameshift SNP in NOD2
13	20763686	20763686	G	-	comments: rs1801002 (del35G), a frameshift mutation in GJB2, associated with hearing loss
13	20797176	21105944	0	-	comments: a 342kb deletion encompassing GJB6, associated with hearing loss

```

Reference files are pre-installed in **$ANNOVAR\_DATA/{build}**, where **{build}** can be either
 hg18, hg19 or hg38. If other builds are needed, contact staff@hpc.nih.gov. To list all builds currently available, type




```
ls $ANNOVAR_DATA
```

Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive -t 4 --gres=lscratch:20**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load annovar**

[user@cn3144 ~]$ **cp $ANNOVAR\_HOME/example/ex1.avinput .**

[user@cn3144 ~]$ **annotate\_variation.pl --geneanno --dbtype refGene --buildver hg38 ex1.avinput $ANNOVAR\_DATA/hg38**

[user@cn3144 ~]$ **table\_annovar.pl ex1.avinput $ANNOVAR\_DATA/hg38 \
 --tempdir /lscratch/$SLURM\_JOB\_ID \
 --thread $SLURM\_CPUS\_ON\_NODE \
 --buildver hg38 \
 --outfile ex1.out \
 --remove \
 --protocol gene,clinvar\_20220320,cosmic70,ljb26\_all,avsnp150,cadd\_1.5 \
 --operation g,f,f,f,f,f \
 --nastring ''**
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. annovar.sh). For example:



```

#!/bin/bash
set -e
module load annovar
annotate_variation.pl --geneanno --dbtype gene --buildver hg38 ex1.avinput $ANNOVAR_DATA/hg38

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] annovar.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. annovar.swarm). For example:



```

convert2annovar.pl -format vcf4 file1.vcf > file1.inp; annotate_variation.pl --geneanno --dbtype gene --buildver hg38 file1.inp $ANNOVAR_DATA/hg38
convert2annovar.pl -format vcf4 file2.vcf > file2.inp; annotate_variation.pl --geneanno --dbtype gene --buildver hg38 file2.inp $ANNOVAR_DATA/hg38
convert2annovar.pl -format vcf4 file3.vcf > file3.inp; annotate_variation.pl --geneanno --dbtype gene --buildver hg38 file3.inp $ANNOVAR_DATA/hg38
convert2annovar.pl -format vcf4 file4.vcf > file4.inp; annotate_variation.pl --geneanno --dbtype gene --buildver hg38 file4.inp $ANNOVAR_DATA/hg38

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f annovar.swarm [-g #] [-t #] --module annovar
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module annovar Loads the annovar module for each subjob in the swarm 
 | |
 | |
 | |






