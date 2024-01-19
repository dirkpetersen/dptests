

document.querySelector('title').textContent = 'PennCNV: Copy Number Variation (CNV) detection from SNP genotyping arrays ';
**PennCNV: Copy Number Variation (CNV) detection from SNP genotyping arrays** 


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



PennCNV is a free software tool for Copy Number Variation (CNV) detection 
from SNP genotyping arrays. Currently it can handle signal intensity data 
from Illumina and Affymetrix arrays. With appropriate preparation of file format, 
it can also handle other types of SNP arrays and oligonucleotide arrays.



### References:


* Wang K, Li M, Hadley D, Liu R, Glessner J, Grant S, Hakonarson H, Bucan M.   

*PennCNV: an integrated hidden Markov model designed for high-resolution 
 copy number variation detection in whole-genome SNP genotyping data.*    

[Genome Research 17:1665-1674, 2007](https://genome.cshlp.org/content/17/11/1665.short).
* Diskin SJ, Li M, Hou C, Yang S, Glessner J, Hakonarson H, Bucan M, Maris JM, Wang K.   

 *Adjustment of genomic waves in signal intensities from whole-genome SNP genotyping platforms.*    

[Nucleic Acids Research 36:e126, 2008](https://academic.oup.com/nar/article/36/19/e126/2409936?login=true).
* Wang K, Chen Z, Tadesse MG, Glessner J, Grant SFA, Hakonarson H, Bucan M, Li M.   

*Modeling genetic inheritance of copy number variations.*    

[Nucleic Acids Research 36:e138, 2008](https://academic.oup.com/nar/article/36/21/e138/2409932?login=true).


Documentation
* [PennCNV Github page](https://github.com/WGLab/PennCNV)
* [PennCNV Home page](http://penncnv.openbioinformatics.org/en/latest/)


Important Notes
* Module Name: PennCNV (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set
	+ **PENNCNV\_HOME**  installation directory
	+ **PENNCNV\_BIN**  executable directory
	+ **PENNCNV\_SRC**  source code directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive** 
[user@cn0862 ~]$ **module load penncnv**
[+] Loading penncnv  1.0.5    

```


```

[user@cn0862 ~]$ **detect\_cnv.pl -h**
Usage:
     detect_cnv.pl [arguments] 

 Optional arguments:
 -v, --verbose use verbose output
 -h, --help print help message
 -m, --man print complete documentation

 Analysis Type:
 --train train optimized HMM model (not recommended to use)
 --test test HMM model to identify CNV
 --wgs test HMM model to identify CNV from transformed wgs data
 --trio posterior CNV calls for father-mother-offspring trio
 --quartet posterior CNV calls for quartet
 --joint joint CNV calls for trio
 --cctest case-control comparison of per-marker CNV frequency
 --validate validate copy number at a pre-specified region

 Input/output Files:
 --listfile <file> a list file containing path to files to be processed
 --output <file> specify output root filename
 --logfile <file> write notification/warning messages to this file
 --hmmfile <file> HMM model file
 --pfbfile <file> population frequency for B allelel file
 --sexfile <file> a 2-column file containing filename and sex (male/female) for chrx CNV calling
 --cnvfile <file> specify CNV call file for use in family-based CNV calling by -trio or -quartet
 --directory <string> specify the directory where signal files are located
 --refchr <string> specify a chromosome for wave adjustment (default: 11 for human)
 --refgcfile <file> a file containing GC percentage of each 1M region of the chromosome specified by --refchr (default: 11 for human)
 --gcmodelfile <file> a file containing GC model for wave adjustment
 --phenofile <file> a file containing phenotype information for each input file for -cctest operation

 CNV output control:
 --minsnp <int> minimum number of SNPs within CNV (default=3)
 --minlength <int> minimum length of bp within CNV
 --minconf <float> minimum confidence score of CNV
 --confidence calculate confidence for each CNV
 --chrx use chrX-specific treatment
 --chry use chrY-specific treatment (available soon)

 Validation-calling arguments:
 --startsnp <string> start SNP of a pre-specified region for --validate operation
 --endsnp <string> end SNP of a pre-specified region for --validate operation
 --delfreq <float> prior deletion frequency of a pre-specified region for --validate operation
 --dupfreq <float> prior duplication frequency of a pre-specified region for --validate operation
 --backfreq <float> background CNV probability for any loci (default: 0.0001)
 --candlist <file> a file containing all candidate CNV regions to be validated

 Misc options
 --loh detect copy-neutral LOH (obselete argument; for SNP arrays only!)
 --exclude\_heterosomic empirically exclude CNVs in heterosomic chromosomes
 --fmprior <numbers> prior belief on CN state for regions with CNV calls
 --denovo\_rate <float> prior belief on genome-wide de novo event rate (default=0.0001)
 --tabout use tab-delimited output
 --coordinate\_from\_input get marker coordindate information from input signal file
 --control\_label <string> the phenotype label for control subjects in the phenotype file (default=control)
 --onesided performed one-sided test for --cctest operation
 --type\_filter <dup|del> used together with --cctest to specify types of CNVs to be tested
 --(no)medianadjust adjust genome-wide LRR such that median=0 (default=ON)
 --(no)bafadjust adjust genome-wide BAF such that median=0.5 (default=ON)
 --(no)sdadjust adjust SD of hidden Markov model based on input signal (default=ON)
 --(no)flush flush input/output buffer (default=ON)
 --bafxhet <float> minimum BAF het rate to predict female gender when -sexfile is not supplied (default=0.1)

 Function: generate CNV calls from high-density SNP genotyping data that
 contains Log R Ratio and B Allele Frequency for each SNP or CN marker. Use -m
 argument to read the complete manual.

 Example: detect\_cnv.pl -test -hmm hhall.hmm -pfb hhall.pfb file1 file2 file3 -log logfile -out outfile
 detect\_cnv.pl -trio -hmm hhall.hmm -pfb hhall.pfb -cnv outfile file1 file2 file3
 detect\_cnv.pl -validate -hmm hhall.hmm -pfb hhall.pfb -startsnp rs100 -endsnp rs200 -delfreq 0.2 file1 file2 file3

 Version: $LastChangedDate: 2013-02-08 11:10:50 -0800 (Fri, 08 Feb 2013) 
...
[user@cn0862 ~]$ **convert\_cnv.pl -h**
...
[user@cn0862 ~]$ **compare\_cnv.pl -h**
...
[user@cn0862 ~]$ **cal\_gc\_snp.pl -h** 
...
[user@cn0862 ~]$ **compile\_pfb.pl -h**
...
[user@cn0862 ~]$ **filter\_cnv.pl -h**
...
[user@cn0862 ~]$ **genomic\_wave.pl -h**
...
[user@cn0862 ~]$ **infer\_snp\_allele.pl -h**
...
[user@cn0862 ~]$ **scan\_region.pl -h**
[user@cn0862 ~]$ **visualize\_cnv.pl -h**
...

```

Exit the application:   


```

[user@cn0862 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```





