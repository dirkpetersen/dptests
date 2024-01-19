

document.querySelector('title').textContent = 'VEP on Biowulf';
VEP on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Plugins](#plugins) 
[Interactive job](#int) 
[Batch job](#sbatch) 
[Swarm of jobs](#swarm) 
 |



VEP (Variant Effect Predictor) determines the effect of your variants (SNPs, insertions, deletions, CNVs or structural variants) on genes, transcripts, and protein sequence, as well as regulatory regions.



### References:


* [McLaren W, Gil L, Hunt SE, Riat HS, Ritchie GR, Thormann A, Flicek P, Cunningham F. The Ensembl Variant Effect Predictor. *Genome Biol.* 2016 Jun 6;17(1):122.](https://www.ncbi.nlm.nih.gov/pubmed/27268795)


Documentation
* [VEP Main Page](http://useast.ensembl.org/info/docs/tools/vep/script/index.html)
* [VEP Tutorial](http://useast.ensembl.org/info/docs/tools/vep/script/vep_tutorial.html)
* [ensembl-vep Git Page](https://github.com/Ensembl/ensembl-vep)
* + [Old VEP versions](#)
		- [VEP Main Page (version 110)](http://jul2023.archive.ensembl.org/info/docs/tools/vep/script/index.html)
		- [VEP Main Page (version 109)](http://feb2023.archive.ensembl.org/info/docs/tools/vep/script/index.html)
		- [VEP Main Page (version 108)](http://oct2022.archive.ensembl.org/info/docs/tools/vep/script/index.html)
		- [VEP Main Page (version 107)](http://jul2022.archive.ensembl.org/info/docs/tools/vep/script/index.html)
		- [VEP Main Page (version 106)](http://apr2022.archive.ensembl.org/info/docs/tools/vep/script/index.html)
		- [VEP Main Page (version 105)](http://dec2021.archive.ensembl.org/info/docs/tools/vep/script/index.html)
		- [VEP Main Page (version 104)](http://may2021.archive.ensembl.org/info/docs/tools/vep/script/index.html)
		- [VEP Main Page (version 103)](http://feb2021.archive.ensembl.org/info/docs/tools/vep/script/index.html)
		- [VEP Main Page (version 102)](http://nov2020.archive.ensembl.org/info/docs/tools/vep/script/index.html)
		- [VEP Main Page (version 101)](http://aug2020.archive.ensembl.org/info/docs/tools/vep/script/index.html)
		- [VEP Main Page (version 100)](http://apr2020.archive.ensembl.org/info/docs/tools/vep/script/index.html)


Important Notes
* Module Name: VEP (see [the modules page](/apps/modules.html) for more information)
* Singlethreaded
* Environment variables set 
	+ VEP\_HOME
	+ VEP\_EXAMPLES
	+ VEP\_CACHEDIR -- reference data and plugins
	+ GS -- location of the GeneSplicer application
	+ PERL5LIB* Example files in $VEP\_EXAMPLES* Reference data in $VEP\_CACHEDIR



By default VEP requires internet connectivity to the Ensembl databases. **THIS IS NOT POSSIBLE ON THE BIOWULF CLUSTER!** Instead, the databases have been locally cached into a version-specific directory **$VEP\_CACHEDIR**, as set by the VEP module, allowing for offline analysis.




Without supplying **--dir $VEP\_CACHEDIR**, VEP will attempt to create a cache directory in your home directory as ~/.vep. **This may cause your /home directory to fill up**, causing unexpected failures to your jobs. Make sure to include **--cache --dir $VEP\_CACHEDIR** with all VEP commands.


Older versions of VEP used **--dir\_cache $VEP\_CACHEDIR** rather than **--dir $VEP\_CACHEDIR**. Both are valid.


### Commands


* convert\_cache.pl
* filter\_vep: filter the output of VEP
* [vep](https://github.com/Ensembl/ensembl-vep#vep): predict effect of variants
* [haplo](https://github.com/Ensembl/ensembl-vep#haplo): Ensembl transcript haplotypes view
* [variant\_recoder](https://github.com/Ensembl/ensembl-vep#recoder): translating between different variant encodings


### Reference Files


**--assembly** is needed for human sequences because there are two available (GRCh37 and GRCh38). For cat, dog, and mouse, no assembly required.


* **--fasta $VEP\_CACHEDIR/cat.fa --species cat**
* **--fasta $VEP\_CACHEDIR/dog.fa --species dog**
* **--fasta $VEP\_CACHEDIR/mouse.fa --species mouse**
* **--fasta $VEP\_CACHEDIR/zebrafish.fa --species zebrafish**
* **--fasta $VEP\_CACHEDIR/GRCh38.fa --species human --assembly GRCh38**
* **--fasta $VEP\_CACHEDIR/GRCh37.fa --species human --assembly GRCh37**


### Plugins


There are a large number of plugins available for use with VEP. Some of these plugins require third-party reference data. Most of this data is available within $VEP\_CACHEDIR, but some are available in
the /fdb tree. Here is an example for using plugins:



```

module load VEP/109
vep \
 -i ${VEP_EXAMPLES}/homo_sapiens_GRCh38.vcf \
 -o example.out \
 --offline \
 --cache \
 --force_overwrite \
 --dir ${VEP_CACHEDIR} \
 --species human \
 --assembly GRCh38 \
 --fasta ${VEP_CACHEDIR}/GRCh38.fa \
 --everything \
 --plugin CSN \
 --plugin Blosum62 \
 --plugin Carol \
 --plugin Condel,$VEP_CACHEDIR/Plugins/config/Condel/config,b \
 --plugin Phenotypes,file=${VEP_CACHEDIR}/Plugins/Phenotypes.pm_human_109_GRCh38.gvf.gz \
 --plugin ExAC,$VEP_CACHEDIR/ExAC.r0.3.sites.vep.vcf.gz \
 --plugin GeneSplicer,$GS/bin/genesplicer,$GS/human,context=200 \
 --plugin CADD,${VEP_CACHEDIR}/CADD_1.4_GRCh38_whole_genome_SNVs.tsv.gz,${VEP_CACHEDIR}/CADD_1.4_GRCh38_InDels.tsv.gz \
 --plugin Downstream \
 --plugin LoFtool \
 --plugin FATHMM,"python ${VEP_CACHEDIR}/fathmm.py"

```

Here is a synopsis of --plugin examples using GRCh38:


* **AlphaMissense**:  


```
--plugin AlphaMissense,file=${VEP_CACHEDIR}/AlphaMissense_hg38.tsv.gz,cols=all
```
* **Blosum62**:  


```
--plugin Blosum62
```
* **CADD**:  


```
--plugin CADD,${VEP_CACHEDIR}/CADD_1.4_GRCh38_whole_genome_SNVs.tsv.gz,${VEP_CACHEDIR}/CADD_1.4_GRCh38_InDels.tsv.gz
```
* **Carol**:  


```
--plugin Carol
```
* **Condel**:  


```
--plugin Condel,${VEP_CACHEDIR}/Plugins/Condel/config,b
```
* **CSN**  


```
--plugin CSN
```
* **dbNSFP**:  


```
--plugin dbNSFP,${VEP_CACHEDIR}/dbNSFP4.2a_hg38.txt.gz,LRT_score,GERP++_RS
```
* **dbscSNV**:  


```
--plugin dbscSNV,${VEP_CACHEDIR}/dbscSNV1.1.txt.gz
```
* **DisGeNET**:  


```
--plugin DisGeNET,file=${VEP_CACHEDIR}/all_variant_disease_pmid_associations_sorted.tsv.gz,disease=1,filter_source=GWASDB&GWASCAT
```
* **Downstream**:  


```
--plugin Downstream
```
* **ExAC**:  


```
--plugin ExAC,${VEP_CACHEDIR}/ExAC.r0.3.1.sites.vep.vcf.gz
```
* **FATHMM**:  


```
--plugin FATHMM,"python ${VEP_CACHEDIR}/fathmm.py"
```
* **G2P**:  


```
--plugin G2P,file=${VEP_CACHEDIR}/SkinG2P.csv
```
* **GeneSplicer**:  


```
--plugin GeneSplicer,$GS/bin/genesplicer,$GS/human,context=200
```
* **HGVSReferenceBase**:  


```
--hgvs --plugin HGVSReferenceBase
```
* **Lof**:  


```
--plugin LoF,loftee_path:${VEP_CACHEDIR}/Plugins/loftee_GRCh38, \
human_ancestor_fa:${VEP_CACHEDIR}/Plugins/loftee_GRCh38/human_ancestor.fa.gz, \
conservation_file:${VEP_CACHEDIR}/Plugins/loftee_GRCh38/loftee.sql, \
gerp_bigwig:${VEP_CACHEDIR}/Plugins/loftee_GRCh38/gerp_conservation_scores.homo_sapiens.GRCh38.bw
```
* **LoFtool**:  


```
--plugin LoFtool
```
* **MTR**:  


```
--plugin MTR,${VEP_CACHEDIR}/mtrflatfile_2.0.txt.gz
```
* **Phenotypes**:  


```
--plugin Phenotypes,file=${VEP_CACHEDIR}/Plugins/Phenotypes.pm_human_109_GRCh38.gvf.gz
```
* **REVEL**:  


```
--plugin REVEL,${VEP_CACHEDIR}/revel_GRCh38_1.3.tsv.gz
```
* **SpliceAI**:  


```
--plugin SpliceAI,snv=${VEP_CACHEDIR}/spliceai_scores.masked.snv.hg38.vcf.gz,indel=${VEP_CACHEDIR}/spliceai_scores.masked.indel.hg38.vcf.gz
```
* **StructuralVariantOverlap**:  


```
--plugin StructuralVariantOverlap,file=${VEP_CACHEDIR}/gnomad_v2_sv.sites.vcf.gz
```


Some of the plugins have multiple versions of reference files.


**NOTE:** The name of the variation file for the Phenotypes plugin depends on the version of VEP being used. In the above example, VEP/**109** is being used, and the variation file is ${VEP\_CACHEDIR}/Plugins/Phenotypes.pm\_homo\_sapiens\_**109**\_GRCh38.gvf.gz. Please run the command **ls -l ${VEP\_CACHEDIR}** to see all the reference files available.


For more information about plugins, type



```
perldoc $VEP_CACHEDIR/Plugins/**[name]**.pm
```

where **[name]** is the name of the plugin.


### LOFTEE Plugin


The [LOFTEE (Loss-Of-Function Transcript Effect Estimator) plugin](https://github.com/konradjk/loftee) is a bit trickier than others, and is not well documented. It is available for VEP versions >= 101. Here are examples for annotating GRCh37 and GRCh38 aligned VCF files:



```

**ml VEP/101**

export PERL5LIB=$PERL5LIB:${VEP_CACHEDIR}/Plugins/loftee_GRCh37
vep \
  --offline --cache --dir ${VEP_CACHEDIR} \
  --input_file ${VEP_EXAMPLES}/homo_sapiens_GRCh37.vcf \
  --species human --assembly GRCh37 --fasta ${VEP_CACHEDIR}/GRCh37.fa \
  --output_file GRCh37.txt --stats_file GRCh37_summary.txt \
  --plugin LoF,loftee_path:${VEP_CACHEDIR}/Plugins/loftee_GRCh37,\
conservation_file:${VEP_CACHEDIR}/Plugins/loftee_GRCh37/phylocsf_gerp.sql,\
human_ancestor_fa:${VEP_CACHEDIR}/Plugins/loftee_GRCh37/human_ancestor.fa.gz >& GRCh37.out

export PERL5LIB=$PERL5LIB:${VEP_CACHEDIR}/Plugins/loftee_GRCh38
vep \
  --offline --cache --dir ${VEP_CACHEDIR} \
  --input_file ${VEP_EXAMPLES}/homo_sapiens_GRCh38.vcf \
  --species human --assembly GRCh38 --fasta ${VEP_CACHEDIR}/GRCh38.fa \
  --output_file GRCh38.txt --stats_file GRCh38_summary.txt \
  --plugin LoF,loftee_path:${VEP_CACHEDIR}/Plugins/loftee_GRCh38,\
human_ancestor_fa:${VEP_CACHEDIR}/Plugins/loftee_GRCh38/human_ancestor.fa.gz,\
conservation_file:${VEP_CACHEDIR}/Plugins/loftee_GRCh38/loftee.sql,\
gerp_bigwig:${VEP_CACHEDIR}/Plugins/loftee_GRCh38/gerp_conservation_scores.homo_sapiens.GRCh38.bw >& GRCh38.out

```

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

[user@cn3144 ~]$ module load VEP
[user@cn3144 ~]$ ln -s $VEP_EXAMPLES/homo_sapiens_GRCh38.vcf .
[user@cn3144 ~]$ vep -i homo_sapiens_GRCh38.vcf -o test.out --offline --cache --dir $VEP_CACHEDIR --species human --assembly GRCh38 --fasta $VEP_CACHEDIR/GRCh38.fa 

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. VEP.sh). For example:



```

#!/bin/bash
module load VEP
vep -i trial1.vcf --offline --cache --dir $VEP_CACHEDIR --fasta $VEP_CACHEDIR/GRCh38.fa --output trial1.out

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] VEP.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. VEP.swarm). For example:



```

vep -i trial1.vcf --offline --cache --dir $VEP_CACHEDIR --fasta $VEP_CACHEDIR/GRCh38.fa --output trial1.out
vep -i trial2.vcf --offline --cache --dir $VEP_CACHEDIR --fasta $VEP_CACHEDIR/GRCh38.fa --output trial2.out
vep -i trial3.vcf --offline --cache --dir $VEP_CACHEDIR --fasta $VEP_CACHEDIR/GRCh38.fa --output trial3.out
vep -i trial4.vcf --offline --cache --dir $VEP_CACHEDIR --fasta $VEP_CACHEDIR/GRCh38.fa --output trial4.out

```

By default, vep will write to the same output file ("variant\_effect\_output.txt") unless directed to do otherwise using the --output option. For swarms of multiple runs, be sure to include this option.


Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f VEP.swarm [-g #] [-t #] --module VEP
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module VEP Loads the VEP module for each subjob in the swarm 
 | |
 | |
 | |






