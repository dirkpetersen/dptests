

document.querySelector('title').textContent = 'CADD: scoring the deleteriousness of SNPs and indels in the Human Genome.';
CADD: scoring the deleteriousness of SNPs and indels in the Human Genome.


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



CADD (Combined Annotation Dependent Depletion) is a tool for scoring the deleteriousness of single nucleotide variants as well as insertion/deletions variants in the human genome. Currently, it supports the builds: GRCh37/hg19 and GRCh38/hg38.



### References:


* Martin Kircher, Daniela M Witten, Preti Jain, Brian J O'Roak, Gregory M Cooper & Jay Shendure   

 *A general framework for estimating the relative pathogenicity of human genetic variants.*  

[Nature Genetics](https://www.nature.com/articles/ng.2892?message-global=remove&page=18)  2014 Feb 2. doi: 10.1038/ng.2892.
* Rentzsch P, Witten D, Cooper GM, Shendure J, Kircher M.   

 *CADD: predicting the deleteriousness of variants throughout the human genome.*  

[Nucleic Acids Res.](https://academic.oup.com/nar/article/47/D1/D886/5146191?login=true)  2018 Oct 29. doi: 10.1093/nar/gky1016.


Documentation
* [CADD Github page](https://github.com/kircherlab/CADD-scripts)


Important Notes
* Module Name: CADD (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set
	+ **CADD\_HOME**  CADD installation directory
	+ **CADD\_BIN**       CADD executable directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive** 
[user@cn4206 ~]$ **module load cadd** 
[+] Loading snakemake  6.5.3
[+] Loading CADD 1.6.post1 on cn4206
[user@cn4206 ~]$ **which CADD**
/usr/local/apps/CADD/1.6.post1/bin/CADD
user@cn4206 ~]$ **CADD -h**
CADD [-o ] [-g ] [-v ] [-a]  -- CADD version 1.6

where:
 -h show this help text
 -o out tsv.gz file (generated from input file name if not set)
 -g genome build (supported are GRCh37 and GRCh38 [default: GRCh38])
 -v CADD version (only v1.6 possible with this set of scripts [default: v1.6])
 -a include annotation in output
 input vcf of vcf.gz file (required)
 -q print basic information about snakemake run
 -p print full information about the snakemake run
 -c number of cores that snakemake is allowed to use [default: 1]
[user@cn4206 ~]$ **CADD -o my\_out.tsv.gz -g GRCh38 $CADD\_TEST/input.vcf**
CADD-v1.6 (c) University of Washington, Hudson-Alpha Institute for Biotechnology and Berlin Institute of Health 2013-2020. All rights reserved.
OUTFILE=my\_out.tsv.gz
Running snakemake pipeline:
snakemake /tmp/tmp.jfNNBZaoNk/input.tsv.gz --use-conda --conda-prefix /vf/db/CADD/1.6.post1/envs --cores 1
--configfile /vf/db/CADD/1.6.post1/config/config\_GRCh38\_v1.6\_noanno.yml --snakefile /vf/db/CADD/1.6.post1/Snakefile -q
Job stats:
job count min threads max threads
---------- ------- ------------- -------------
annotation 1 1 1
imputation 1 1 1
join 1 1 1
prepare 1 1 1
prescore 1 1 1
score 1 1 1
total 6 1 1

Opening /vf/db/CADD/1.6.post1/data/prescored/GRCh38\_v1.6/no\_anno/gnomad.genomes.r3.0.indel.tsv.gz...
Opening /vf/db/CADD/1.6.post1/data/prescored/GRCh38\_v1.6/no\_anno/whole\_genome\_SNVs.tsv.gz...
Possible precedence issue with control flow operator at /vf/db/CADD/1.6.post1/envs/e507b9d6fe0a5d0ae12148a5d819bc0e/lib/site\_perl/5.26.2/Bio/DB/IndexedBase.pm line 805.

CADD scored variants written to file: my\_out.tsv.gz

[user@cn4206 ~]$ **CADD -o my\_out.tsv.gz -g GRCh38 $CADD\_TEST/input.vcf**
CADD-v1.6 (c) University of Washington, Hudson-Alpha Institute for Biotechnology and Berlin Institute of Health 2013-2020. All rights reserved.
Running snakemake pipeline:
snakemake /tmp/tmp.scK9I9hCpe/input.tsv.gz --use-conda --conda-prefix /vf/db/CADD/1.6.post1/envs --cores 1
--configfile /vf/db/CADD/1.6.post1/config/config\_GRCh37\_v1.6\_noanno.yml --snakefile /vf/db/CADD/1.6.post1/Snakefile -q
Job stats:
job count min threads max threads
---------- ------- ------------- -------------
annotation 1 1 1
imputation 1 1 1
join 1 1 1
prepare 1 1 1
prescore 1 1 1
score 1 1 1
total 6 1 1

Opening /vf/db/CADD/1.6.post1/data/prescored/GRCh37\_v1.6/no\_anno/InDels.tsv.gz...
Opening /vf/db/CADD/1.6.post1/data/prescored/GRCh37\_v1.6/no\_anno/whole\_genome\_SNVs.tsv.gz...
Possible precedence issue with control flow operator at /vf/db/CADD/1.6.post1/envs/e507b9d6fe0a5d0ae12148a5d819bc0e/lib/site\_perl/5.26.2/Bio/DB/IndexedBase.pm line 805.

CADD scored variants written to file: my\_out.tsv.gz

```

End the interactive session:

```

[user@cn4206 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```





