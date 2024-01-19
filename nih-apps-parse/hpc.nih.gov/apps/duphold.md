

document.querySelector('title').textContent = "Duphold";
DUPHOLD on Biowulf


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



Duphold: scalable, depth-based annotation and curation of high-confidence structural variant calls.



### References:


* Pedersen BS, Quinlan AR.
 [**Duphold: scalable, depth-based annotation and curation of high-confidence structural variant calls.**](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6479422/)
*Gigascience. 2019 Apr 1;8(4):giz040. doi: 10.1093/gigascience/giz040. PMID: 31222198; PMCID: PMC6479422.*


Documentation
* [DUPHOLD Main Site](https://github.com/brentp/duphold)


Important Notes
* Module Name: duphold (see [the modules page](/apps/modules.html) for more information)
 * Multithreaded
* Environment variables set 
	+ DUPHOLD\_HOME



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive**
[user@cn4338 ~]$ **module load duphold**
[+] Loading duphold  0.2.3  on cn4338 

```


Example
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Running help command:



```

[user@cn4338 data]$ **duphold --help**
version: 0.2.3

  Usage: duphold [options]

Options:
  -v --vcf  path to sorted SV VCF/BCF
 -b --bam  path to indexed BAM/CRAM
 -f --fasta  indexed fasta reference.
 -s --snp  optional path to snp/indel VCF/BCF with which to annotate SVs. BCF is highly recommended as it's much faster to parse.
 -t --threads  number of decompression threads. [default: 4]
 -o --output  output VCF/BCF (default is VCF to stdout) [default: -]
 -d --drop drop all samples from a multi-sample --vcf \*except\* the sample in --bam. useful for parallelization by sample followed by merge.
 -h --help show help

```

Annotate an SV:



```

[user@cn4338] **cp -a /usr/local/apps/duphold/0.2.3/test\_data .**
[user@cn4338 test_data]$ **duphold \
 --threads 4 \
 --vcf sparse\_in.vcf \
 --bam sparse.cram \
 --fasta sparse.fa \
 --output output.bcf** 
#To view output, load samtools and view with bcftools
[user@cn4338 test_data] **module load samtools**
[user@cn4338 test_data] **bcftools view test-out.bcf**
##fileformat=VCFv4.2
...
##bcftools_viewVersion=1.4-19-g1802ff3+htslib-1.4-29-g42bfe70
##bcftools_viewCommand=view CHM1_CHM13/full.37d5.vcf.gz; Date=Mon Sep 24 13:48:04 2018
...
##bcftools_viewVersion=1.17+htslib-1.17
##bcftools_viewCommand=view test-out.bcf; Date=Thu May 25 12:49:34 2023
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Eluc-CR2.F
NW_017858824.1  135118  72454   N       DEL   5875.46 .       SVTYPE=DEL;END=135332;CIPOS=0,0;CIEND=0,0;CIPOS95=0,0;CIEND95=0,0;GCF=0.306977  GT:DP:DHFC:DHFFC:DHBFC:DHSP     0/1:200:1.91667:0.597403:1.76923:0

```



|  |  |
| --- | --- |
|  For more information on pre and post processing, please visit the [Duphold Github Page](https://github.com/brentp/duphold) | |








