

document.querySelector('title').textContent = 'Clinical interpretation of genetic variants by the 2015 ACMG-AMP guideline on Biowulf';
**Clinical interpretation of genetic variants by the 2015 ACMG-AMP guideline on Biowulf**


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



InterVar is bioinformatics software tool 
for clinical interpretation of genetic variants by the ACMG-AMP 2015 guidelines,
which are published by the American College of Medical Genetics and Genomics (ACMG) 
and the Association for Molecular Pathology (AMP).
These sre the standards for the clinical interpretation of sequence variants 
with respect to human diseases on the basis of 28 criteria. 



### References:


* Quan Li and Kai Wang,
InterVar: Clinical Interpretation of Genetic Variants by the 2015 ACMG-AMP Guidelines,
The American Journal of Human Genetics, 2017, Volume 100, Issue 2, Pages 267-280


Documentation
* InterVar Main Site: [InterVar on github](https://github.com/WGLab/InterVar)


Important Notes
* Module Name: InterVar (see [the modules page](/apps/modules.html) for more information)
* Singlethreaded
* Unusual environment variables set
	+ **INTERVAR\_HOME**   InterVar installation directory
	+ **INTERVAR\_BIN**         InterVar executable directory
	+ **INTERVAR\_DATA**     InterVar data directory
	+ **INTERVAR\_TEST**      InterVar examples directory
	+ **ANNOVAR\_HOME**     ANNOVAR installation directory
	+ **ANNOVAR\_DATA**       ANNOVAR data directory* Example files in **$INTERVAR\_TEST**



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=16g**
salloc.exe: Pending job allocation 59748321
salloc.exe: job 59748321 queued and waiting for resources
salloc.exe: job 59748321 has been allocated resources
salloc.exe: Granted job allocation 59748321
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3624 are ready for job
[user@cn3144 ~]$ **module load InterVar**
[+] Loading InterVar 2.1.3 ...
[user@cn3144 ~]$ **InterVar -i $INTERVAR\_TEST/ex1.avinput -d $ANNOVAR\_DATA/hg19 -o result\_hg19**
=============================================================================
InterVar
Interpretation of Pathogenic/Benign for variants using python scripts.

.####.##....##.########.########.########..##.....##....###....########.
..##..###...##....##....##.......##.....##.##.....##...##.##...##.....##
..##..####..##....##....##.......##.....##.##.....##..##...##..##.....##
..##..##.##.##....##....######...########..##.....##.##.....##.########.
..##..##..####....##....##.......##...##....##...##..#########.##...##..
..##..##...###....##....##.......##....##....##.##...##.....##.##....##.
.####.##....##....##....########.##.....##....###....##.....##.##.....##

=============================================================================

%prog 2.0.2 20190327
Written by Quan LI,leequan@gmail.com.
InterVar is free for non-commercial use without warranty.
Please contact the authors for commercial use.
Copyright (C) 2016 Wang Genomic Lab
============================================================================

Notice: Your command of InterVar is ['/usr/local/apps/InterVar/2.1.3/bin/Intervar.py', '-t', '/usr/local/apps/InterVar/DOWNLOADS/InterVar/intervardb', '--annotate_variation=/usr/local/apps/annovar/2018-04-16/annotate_variation.pl', '--table_annovar=/usr/local/apps/annovar/2018-04-16/table_annovar.pl', '--convert2annovar=/usr/local/apps/annovar/2018-04-16/convert2annovar.pl', '-i', '/usr/local/apps/InterVar/TEST_DATA/example/ex1.avinput', '-d', '/fdb/annovar/2018-04-16/hg19', '-o', 'result_hg19']
INFO: The options are {'pp2_genes': '/usr/local/apps/InterVar/DOWNLOADS/InterVar/intervardb/PP2.genes.hg19', 'inputfile': '/usr/local/apps/InterVar/TEST_DATA/example/ex1.avinput', 'exclude_snps': '/usr/local/apps/InterVar/DOWNLOADS/InterVar/intervardb/ext.variants.hg19', 'annotate_variation': '/usr/local/apps/annovar/2018-04-16/annotate_variation.pl', 'skip_annovar': False, 'ps4_snps': '/usr/local/apps/InterVar/DOWNLOADS/InterVar/intervardb/PS4.variants.hg19', 'mim_domin': '/usr/local/apps/InterVar/DOWNLOADS/InterVar/intervardb/mim_domin.txt', 'current_version': 'Intervar_20190327', 'bs2_snps': '/usr/local/apps/InterVar/DOWNLOADS/InterVar/intervardb/BS2_hom_het.hg19', 'evidence_file': 'None', 'public_dev': 'https://github.com/WGLab/InterVar/releases', 'otherinfo': 'TRUE', 'database_names': 'refGene esp6500siv2_all 1000g2015aug avsnp147 dbnsfp33a clinvar_20190305 gnomad_genome dbscsnv11 dbnsfp31a_interpro rmsk ensGene knownGene', 'mim_pheno': '/usr/local/apps/InterVar/DOWNLOADS/InterVar/intervardb/mim_pheno.txt', 'table_annovar': '/usr/local/apps/annovar/2018-04-16/table_annovar.pl', 'buildver': 'hg19', 'inputfile_type': 'AVinput', 'onetranscript': 'FALSE', 'mim2gene': '/usr/local/apps/InterVar/DOWNLOADS/InterVar/intervardb/mim2gene.txt', 'orpha': '/usr/local/apps/InterVar/DOWNLOADS/InterVar/intervardb/orpha.txt', 'ps1_aa': '/usr/local/apps/InterVar/DOWNLOADS/InterVar/intervardb/PS1.AA.change.patho.hg19', 'mim_adultonset': '/usr/local/apps/InterVar/DOWNLOADS/InterVar/intervardb/mim_adultonset.txt', 'knowngenecanonical': '/usr/local/apps/InterVar/DOWNLOADS/InterVar/intervardb/knownGeneCanonical.txt.hg19', 'outfile': 'result_hg19', 'morbidmap': '/usr/local/apps/InterVar/DOWNLOADS/InterVar/intervardb/morbidmap', 'convert2annovar': '/usr/local/apps/annovar/2018-04-16/convert2annovar.pl', 'database_locat': '/fdb/annovar/2018-04-16/hg19', 'database_intervar': '/usr/local/apps/InterVar/DOWNLOADS/InterVar/intervardb', 'lof_genes': '/usr/local/apps/InterVar/DOWNLOADS/InterVar/intervardb/PVS1.LOF.genes.hg19', 'disorder_cutoff': '0.01', 'mim_recessive': '/usr/local/apps/InterVar/DOWNLOADS/InterVar/intervardb/mim_recessive.txt', 'pm1_domain': '/usr/local/apps/InterVar/DOWNLOADS/InterVar/intervardb/PM1_domains_with_benigns.hg19', 'mim_orpha': '/usr/local/apps/InterVar/DOWNLOADS/InterVar/intervardb/mim_orpha.txt', 'bp1_genes': '/usr/local/apps/InterVar/DOWNLOADS/InterVar/intervardb/BP1.genes.hg19'}
Warning: the folder of /fdb/annovar/2018-04-16/hg19 is already created!
perl /usr/local/apps/annovar/2018-04-16/table_annovar.pl /usr/local/apps/InterVar/TEST_DATA/example/ex1.avinput /fdb/annovar/2018-04-16/hg19 -buildver hg19 -remove -out result_hg19 -protocol refGene,esp6500siv2_all,1000g2015aug_all,avsnp147,dbnsfp33a,clinvar_20190305,gnomad_genome,dbscsnv11,dbnsfp31a_interpro,rmsk,ensGene,knownGene  -operation  g,f,f,f,f,f,f,f,f,r,g,g   -nastring . --otherinfo
-----------------------------------------------------------------
NOTICE: Processing operation=g protocol=refGene

NOTICE: Running with system command <annotate_variation.pl -geneanno -buildver hg19 -dbtype refGene -outfile result_hg19.refGene -exonsort /usr/local/apps/InterVar/TEST_DATA/example/ex1.avinput /fdb/annovar/2018-04-16/hg19>
NOTICE: Output files were written to result_hg19.refGene.variant_function, result_hg19.refGene.exonic_variant_function
NOTICE: Reading gene annotation from /fdb/annovar/2018-04-16/hg19/hg19_refGene.txt ... Done with 78239 transcripts (including 18578 without coding sequence annotation) for 28293 unique genes
NOTICE: Processing next batch with 18 unique variants in 18 input lines
NOTICE: Reading FASTA sequences from /fdb/annovar/2018-04-16/hg19/hg19_refGeneMrna.fa ... Done with 24 sequences
WARNING: A total of 465 sequences will be ignored due to lack of correct ORF annotation
-----------------------------------------------------------------
NOTICE: Processing operation=f protocol=esp6500siv2_all

NOTICE: Running system command <annotate_variation.pl -filter -dbtype esp6500siv2_all -buildver hg19 -outfile result_hg19 /usr/local/apps/InterVar/TEST_DATA/example/ex1.avinput /fdb/annovar/2018-04-16/hg19>
NOTICE: the --dbtype esp6500siv2_all is assumed to be in generic ANNOVAR database format
NOTICE: Variants matching filtering criteria are written to result_hg19.hg19_esp6500siv2_all_dropped, other variants are written to result_hg19.hg19_esp6500siv2_all_filtered
NOTICE: Processing next batch with 18 unique variants in 18 input lines
NOTICE: Database index loaded. Total number of bins is 594771 and the number of bins to be scanned is 12
NOTICE: Scanning filter database /fdb/annovar/2018-04-16/hg19/hg19_esp6500siv2_all.txt...Done
-----------------------------------------------------------------
NOTICE: Processing operation=f protocol=1000g2015aug_all

NOTICE: Running system command <annotate_variation.pl -filter -dbtype 1000g2015aug_all -buildver hg19 -outfile result_hg19 /usr/local/apps/InterVar/TEST_DATA/example/ex1.avinput /fdb/annovar/2018-04-16/hg19>
NOTICE: Variants matching filtering criteria are written to result_hg19.hg19_ALL.sites.2015_08_dropped, other variants are written to result_hg19.hg19_ALL.sites.2015_08_filtered
NOTICE: Processing next batch with 18 unique variants in 18 input lines
NOTICE: Database index loaded. Total number of bins is 2824642 and the number of bins to be scanned is 17
NOTICE: Scanning filter database /fdb/annovar/2018-04-16/hg19/hg19_ALL.sites.2015_08.txt...Done
-----------------------------------------------------------------
NOTICE: Processing operation=f protocol=avsnp147

NOTICE: Running system command <annotate_variation.pl -filter -dbtype avsnp147 -buildver hg19 -outfile result_hg19 /usr/local/apps/InterVar/TEST_DATA/example/ex1.avinput /fdb/annovar/2018-04-16/hg19>
NOTICE: Variants matching filtering criteria are written to result_hg19.hg19_avsnp147_dropped, other variants are written to result_hg19.hg19_avsnp147_filtered
NOTICE: Processing next batch with 18 unique variants in 18 input lines
NOTICE: Database index loaded. Total number of bins is 27868332 and the number of bins to be scanned is 17
NOTICE: Scanning filter database /fdb/annovar/2018-04-16/hg19/hg19_avsnp147.txt...Done
-----------------------------------------------------------------
NOTICE: Processing operation=f protocol=dbnsfp33a
NOTICE: Finished reading 66 column headers for '-dbtype dbnsfp33a'

NOTICE: Running system command <annotate_variation.pl -filter -dbtype dbnsfp33a -buildver hg19 -outfile result_hg19 /usr/local/apps/InterVar/TEST_DATA/example/ex1.avinput /fdb/annovar/2018-04-16/hg19 -otherinfo>
NOTICE: the --dbtype dbnsfp33a is assumed to be in generic ANNOVAR database format
NOTICE: Variants matching filtering criteria are written to result_hg19.hg19_dbnsfp33a_dropped, other variants are written to result_hg19.hg19_dbnsfp33a_filtered
NOTICE: Processing next batch with 18 unique variants in 18 input lines
NOTICE: Database index loaded. Total number of bins is 550512 and the number of bins to be scanned is 11
NOTICE: Scanning filter database /fdb/annovar/2018-04-16/hg19/hg19_dbnsfp33a.txt...Done
-----------------------------------------------------------------
NOTICE: Processing operation=f protocol=clinvar_20190305
NOTICE: Finished reading 5 column headers for '-dbtype clinvar_20190305'

NOTICE: Running system command <annotate_variation.pl -filter -dbtype clinvar_20190305 -buildver hg19 -outfile result_hg19 /usr/local/apps/InterVar/TEST_DATA/example/ex1.avinput /fdb/annovar/2018-04-16/hg19 -otherinfo>
NOTICE: the --dbtype clinvar_20190305 is assumed to be in generic ANNOVAR database format
NOTICE: Variants matching filtering criteria are written to result_hg19.hg19_clinvar_20190305_dropped, other variants are written to result_hg19.hg19_clinvar_20190305_filtered
NOTICE: Processing next batch with 18 unique variants in 18 input lines
NOTICE: Database index loaded. Total number of bins is 45822 and the number of bins to be scanned is 10
NOTICE: Scanning filter database /fdb/annovar/2018-04-16/hg19/hg19_clinvar_20190305.txt...Done
-----------------------------------------------------------------
NOTICE: Processing operation=f protocol=gnomad_genome
NOTICE: Finished reading 8 column headers for '-dbtype gnomad_genome'

NOTICE: Running system command <annotate_variation.pl -filter -dbtype gnomad_genome -buildver hg19 -outfile result_hg19 /usr/local/apps/InterVar/TEST_DATA/example/ex1.avinput /fdb/annovar/2018-04-16/hg19 -otherinfo>
NOTICE: the --dbtype gnomad_genome is assumed to be in generic ANNOVAR database format
NOTICE: Variants matching filtering criteria are written to result_hg19.hg19_gnomad_genome_dropped, other variants are written to result_hg19.hg19_gnomad_genome_filtered
NOTICE: Processing next batch with 18 unique variants in 18 input lines
NOTICE: Database index loaded. Total number of bins is 28127612 and the number of bins to be scanned is 17
NOTICE: Scanning filter database /fdb/annovar/2018-04-16/hg19/hg19_gnomad_genome.txt...Done
-----------------------------------------------------------------
NOTICE: Processing operation=f protocol=dbscsnv11
NOTICE: Finished reading 2 column headers for '-dbtype dbscsnv11'

NOTICE: Running system command <annotate_variation.pl -filter -dbtype dbscsnv11 -buildver hg19 -outfile result_hg19 /usr/local/apps/InterVar/TEST_DATA/example/ex1.avinput /fdb/annovar/2018-04-16/hg19 -otherinfo>
NOTICE: the --dbtype dbscsnv11 is assumed to be in generic ANNOVAR database format
NOTICE: Variants matching filtering criteria are written to result_hg19.hg19_dbscsnv11_dropped, other variants are written to result_hg19.hg19_dbscsnv11_filtered
NOTICE: Processing next batch with 18 unique variants in 18 input lines
NOTICE: Database index loaded. Total number of bins is 421415 and the number of bins to be scanned is 7
NOTICE: Scanning filter database /fdb/annovar/2018-04-16/hg19/hg19_dbscsnv11.txt...Done
-----------------------------------------------------------------
NOTICE: Processing operation=f protocol=dbnsfp31a_interpro
NOTICE: Finished reading 1 column headers for '-dbtype dbnsfp31a_interpro'

NOTICE: Running system command <annotate_variation.pl -filter -dbtype dbnsfp31a_interpro -buildver hg19 -outfile result_hg19 /usr/local/apps/InterVar/TEST_DATA/example/ex1.avinput /fdb/annovar/2018-04-16/hg19 -otherinfo>
NOTICE: the --dbtype dbnsfp31a_interpro is assumed to be in generic ANNOVAR database format
NOTICE: Variants matching filtering criteria are written to result_hg19.hg19_dbnsfp31a_interpro_dropped, other variants are written to result_hg19.hg19_dbnsfp31a_interpro_filtered
NOTICE: Processing next batch with 18 unique variants in 18 input lines
NOTICE: Database index loaded. Total number of bins is 275326 and the number of bins to be scanned is 5
NOTICE: Scanning filter database /fdb/annovar/2018-04-16/hg19/hg19_dbnsfp31a_interpro.txt...Done
-----------------------------------------------------------------
NOTICE: Processing operation=r protocol=rmsk

NOTICE: Running with system command <annotate_variation.pl -regionanno -dbtype rmsk -buildver hg19 -outfile result_hg19 /usr/local/apps/InterVar/TEST_DATA/example/ex1.avinput /fdb/annovar/2018-04-16/hg19>
NOTICE: Output file is written to result_hg19.hg19_rmsk
NOTICE: Reading annotation database /fdb/annovar/2018-04-16/hg19/hg19_rmsk.txt ... Done with 5298130 regions
NOTICE: Finished region-based annotation on 18 genetic variants
-----------------------------------------------------------------
NOTICE: Processing operation=g protocol=ensGene

NOTICE: Running with system command <annotate_variation.pl -geneanno -buildver hg19 -dbtype ensGene -outfile result_hg19.ensGene -exonsort /usr/local/apps/InterVar/TEST_DATA/example/ex1.avinput /fdb/annovar/2018-04-16/hg19>
NOTICE: Output files were written to result_hg19.ensGene.variant_function, result_hg19.ensGene.exonic_variant_function
NOTICE: Reading gene annotation from /fdb/annovar/2018-04-16/hg19/hg19_ensGene.txt ... Done with 196501 transcripts (including 101155 without coding sequence annotation) for 57905 unique genes
NOTICE: Processing next batch with 18 unique variants in 18 input lines
NOTICE: Reading FASTA sequences from /fdb/annovar/2018-04-16/hg19/hg19_ensGeneMrna.fa ... Done with 22 sequences
WARNING: A total of 6780 sequences will be ignored due to lack of correct ORF annotation
-----------------------------------------------------------------
NOTICE: Processing operation=g protocol=knownGene

NOTICE: Running with system command <annotate_variation.pl -geneanno -buildver hg19 -dbtype knownGene -outfile result_hg19.knownGene -exonsort /usr/local/apps/InterVar/TEST_DATA/example/ex1.avinput /fdb/annovar/2018-04-16/hg19>
NOTICE: Output files were written to result_hg19.knownGene.variant_function, result_hg19.knownGene.exonic_variant_function
NOTICE: Reading gene annotation from /fdb/annovar/2018-04-16/hg19/hg19_knownGene.txt ... Done with 78963 transcripts (including 18502 without coding sequence annotation) for 28495 unique genes
NOTICE: Processing next batch with 18 unique variants in 18 input lines
NOTICE: Reading FASTA sequences from /fdb/annovar/2018-04-16/hg19/hg19_knownGeneMrna.fa ... Done with 47 sequences
WARNING: A total of 43 sequences will be ignored due to lack of correct ORF annotation
-----------------------------------------------------------------
NOTICE: Multianno output file is written to result_hg19.hg19_multianno.txt
Notice: Begin the variants interpretation by InterVar
Notice: About 18 lines in your variant file!
Notice: About 22 variants has been processed by InterVar
Notice: The InterVar is finished, the output file is [ result_hg19.hg19_multianno.txt.intervar ]
=============================================================================
........................................................................
.####.##....##.########.########.########..##.....##....###....########.
..##..###...##....##....##.......##.....##.##.....##...##.##...##.....##
..##..####..##....##....##.......##.....##.##.....##..##...##..##.....##
..##..##.##.##....##....######...########..##.....##.##.....##.########.
..##..##..####....##....##.......##...##....##...##..#########.##...##..
..##..##...###....##....##.......##....##....##.##...##.....##.##....##.
.####.##....##....##....########.##.....##....###....##.....##.##.....##
.......................................................................
Thanks for using InterVar!
Report bugs to leequan@gmail.com;
InterVar homepage: <http://wInterVar.wglab.org>

=============================================================================

[user@cn3144 ~]$ **exit**
exit
salloc.exe: Relinquishing job allocation 59748321
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. InterVar.sh). For example:



```

#!/bin/bash
module load InterVar       
InterVar -i $INTERVAR_TEST/ex1.avinput -d $ANNOVAR_DATA/hg19 -o result_hg19

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] InterVar.sh
```







