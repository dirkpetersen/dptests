

document.querySelector('title').textContent = 'MoChA on Biowulf';
MoChA on Biowulf


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



MoChA is a bcftools extension to call mosaic chromosomal alterations 
starting from phased VCF files 
with either B Allele Frequency (BAF) and Log R Ratio (LRR) or allelic depth (AD).



### References:


* Loh P., Genovese G., McCarroll S., Price A. et al.   

*Insights about clonal expansions from 8,342 mosaic chromosomal alterations.*    

[Nature **559** (2018), 350–355.](https://www.nature.com/articles/s41586-018-0321-x) 
PMID: 29995854; DOI: 10.1038/s41586-018-0321-x
* Loh P., Genovese G. and McCarroll S.   

*Monogenic and polygenic inheritance become instruments for clonal selection*
 [Nature **584** ((2020), 136–141.](https://www.nature.com/articles/s41586-020-2430-6?elqTrackId=57c8cccc30f54cb182de08a725693b0d) PMID: 32581363; DOI: 10.1038/s41586-020-2430-6.


Documentation
* [MoChA Home page](http://software.broadinstitute.org/software/mocha)
* [MoChA Github page](https://github.com/freeseek/mocha)


Important Notes
* Module Name: mocha (see [the modules page](/apps/modules.html) for more information)
* Multithreaded app (set --threads $SLURM\_CPUS\_PER\_TASK)
* Reference data in /fdb/mocha/GRCh37* Unusual environment variables set
	+ **MOCHA\_HOME**  installation directory
	+ **MOCHA\_BIN**  executable directory
	+ **MOCHA\_WDL**  WDL pipeline scripts folder
	+ **MOCHA\_REF**  reference data directory
	+ **MOCHA\_DATA**  sample data directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive -c 8 --mem 20g --gres=lscratch:20 --time=2:00:00**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job


```

[user@cn3144 ~]$ **module load mocha**
+] Loading samtools 1.16  ...
[+] Loading zlib 1.2.11  ...
[+] Loading bedtools  2.29.2
[+] Loading gcc  9.2.0  ...
[+] Loading GSL 2.6 for GCC 9.2.0 ...
[+] Loading openmpi 3.1.4  for GCC 9.2.0
[+] Loading ImageMagick  7.0.8  on cn2385
[+] Loading HDF5  1.10.4
[+] Loading NetCDF 4.7.4_gcc9.2.0
[+] Loading pandoc  2.10.1  on cn2385
[+] Loading pcre2 10.21  ...
[+] Loading R 4.0.0
[+] Loading mocha  1.16

```

Available executables are in stored in the folder pointed to by the environment variable MOCHA_BIN:

```

[user@cn3144 ~]$ **ls $MOCHA\_BIN**
extendFMT  mocha  mochatools  trio-phase

```

To print the usage message for an executable, type its name without arguments. For example:

```

[user@cn3144 ~]$ **mocha**
Genome reference assembly was not specified with --rules or --rules-file

About:   MOsaic CHromosomal Alterations caller, requires phased genotypes (GT)
         and either B-allele frequency (BAF) and Log R Ratio intensity (LRR)
         or allelic depth coverage (AD). (version 2020-09-01 https://github.com/freeseek/mocha)
Usage:   mocha [OPTIONS] <in.vcf>

Required options:
    -r, --rules <assembly>[?]      predefined genome reference rules, 'list' to print available settings, append '?' for details
    -R, --rules-file <file>        genome reference rules, space/tab-delimited CHROM:FROM-TO,TYPE

General Options:
    -x, --sex <file>               file including information about the gender of the samples
        --call-rate <file>         file including information about the call_rate of the samples
    -s, --samples [^]<list>        comma separated list of samples to include (or exclude with "^" prefix)
    -S, --samples-file [^]<file>   file of samples to include (or exclude with "^" prefix)
        --force-samples            only warn about unknown subset samples
    -v, --variants [^]<file>       tabix-indexed [compressed] VCF/BCF file containing variants
    -t, --targets [^]<region>      restrict to comma-separated list of regions. Exclude regions with "^" prefix
    -T, --targets-file [^]<file>   restrict to regions listed in a file. Exclude regions with "^" prefix
    -f, --apply-filters <list>     require at least one of the listed FILTER strings (e.g. "PASS,.")
                                   to include (or exclude with "^" prefix) in the analysis
    -p  --cnp <file>               list of regions to genotype in BED format
        --mhc <region>             MHC region to exclude from analysis (will be retained in the output)
        --kir <region>             KIR region to exclude from analysis (will be retained in the output)
        --threads <int>            number of extra output compression threads [0]

Output Options:
    -o, --output <file>            write output to a file [no output]
    -O, --output-type <b|u|z|v>    b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: uncompressed VCF [v]
        --no-version               do not append version and command line to the header
    -a  --no-annotations           omit Ldev and Bdev FORMAT from output VCF (requires --output)
        --no-log                   suppress progress report on standard error
    -l  --log <file>               write log to file [standard error]
    -m, --mosaic-calls <file>      write mosaic chromosomal alterations to a file [standard output]
    -g, --genome-stats <file>      write sample genome-wide statistics to a file [no output]
    -u, --ucsc-bed <file>          write UCSC bed track to a file [no output]

HMM Options:
        --bdev-LRR-BAF <list>      comma separated list of inverse BAF deviations for LRR+BAF model [-2.0,-4.0,-6.0,10.0,6.0,4.0]
        --bdev-BAF-phase <list>    comma separated list of inverse BAF deviations for BAF+phase model
                                   [6.0,8.0,10.0,15.0,20.0,30.0,50.0,80.0,100.0,150.0,200.0]
        --min-dist <int>           minimum base pair distance between consecutive sites for WGS data [400]
        --adjust-BAF-LRR <int>     minimum number of genotypes for a cluster to median adjust BAF and LRR (-1 for no adjustment) [5]
        --regress-BAF-LRR <int>    minimum number of genotypes for a cluster to regress BAF against LRR (-1 for no regression) [15]
        --LRR-GC-order <int>       order of polynomial to regress LRR against local GC content (-1 for no regression) [2]
        --xy-prob <float>          transition probability [1e-06]
        --err-prob <float>         uniform error probability [1e-02]
        --flip-prob <float>        phase flip probability [1e-02]
        --centromere-loss <float>  penalty to avoid calls spanning centromeres [1e-04]
        --telomere-gain <float>    telomere advantage to prioritize CN-LOHs [1e-02]
        --x-telomere-gain <float>  X telomere advantage to prioritize mLOX [1e-04]
        --y-telomere-gain <float>  Y telomere advantage to prioritize mLOY [1e-05]
        --short-arm-chrs <list>    list of chromosomes with short arms [13,14,15,21,22,chr13,chr14,chr15,chr21,chr22]
        --use-short-arms           use variants in short arms [FALSE]
        --use-centromeres          use variants in centromeres [FALSE]
        --use-no-rules-chrs        use chromosomes without centromere rules  [FALSE]
        --LRR-weight <float>       relative contribution from LRR for LRR+BAF  model [0.2]
        --LRR-hap2dip <float>      difference between LRR for haploid and diploid [0.45]
        --LRR-cutoff <float>       cutoff between LRR for haploid and diploid used to infer gender [estimated from X nonPAR]
...

```

In order to run the mocha executable on sample data, first download the data to your current directory (  **about 3.5GB**  ), then run the sample command below to output results  to the compressed BCF file (  **about 3.5GB,uncompressed VCF file about 14GB**  ):

```

[user@cn3144 ~]$ **cp $MOCHA\_DATA/\* .**
[user@cn3144 ~]$ **mocha -g GRCh37 -p $MOCHA\_REF/cnp.grch37.bed.gz --LRR-GC-order 0 -c calls.tsv -z stats.tsv -u ucsc.bed -Ob --threads 7 -o output.bcf kgp\_om25.8v1-1.as.bcf**
MoChA 2020-09-01 https://github.com/freeseek/mocha
Genome reference: GRCh37
Regions to genotype: /fdb/mocha/GRCh37/cnp.grch37.bed.gz
BAF deviations for LRR+BAF model: -2.0,-4.0,-6.0,10.0,6.0,4.0
BAF deviations for BAF+phase model: 6.0,8.0,10.0,15.0,20.0,30.0,50.0,80.0,100.0,150.0,200.0
Minimum base pair distance between consecutive sites: 400
Order of polynomial in local GC content to be used to regress LRR against GC: 0
Transition probability: 1e-06
Uniform error probability: 0.01
Phase flip probability: 0.01
Centromere penalty: 0.0001
Telomere advantage: 0.01
X telomere advantage: 0.0001
Y telomere advantage: 1e-05
List of short arms: 13,14,15,21,22,chr13,chr14,chr15,chr21,chr22
Use variants in short arms: FALSE
Use variants in centromeres: FALSE
Use chromosomes without centromere rules: FALSE
Relative contribution from LRR for LRR+BAF model: 0.20
Difference between LRR for haploid and diploid: 0.69
Using genome assembly from GRCh37
Loading 2 sample(s) from the VCF file
Read 38 variants from contig 1
Read 18 variants from contig 2
Read 22 variants from contig 3
Read 12 variants from contig 4
Read 10 variants from contig 5
Read 13 variants from contig 6
Read 13 variants from contig 7
Read 7 variants from contig 8
Read 11 variants from contig 9
Read 15 variants from contig 10
...
[user@cn3144 ~]$ **more output.vcf**
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=G,Type=Integer,Description="Allelic Depths of REF and ALT(s) in the order listed">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##contig=<ID=1>
##contig=<ID=2>
##contig=<ID=3>
##contig=<ID=4>
##contig=<ID=5>
##contig=<ID=6>
##contig=<ID=7>
...
##FORMAT=<ID=Bdev_Phase,Number=1,Type=Integer,Description="BAF deviation phase, if available">
##bcftools_pluginVersion=1.10+htslib-1.10
##bcftools_pluginCommand=/usr/local/apps/mocha/1.10.2/lib/mocha.so --LRR-GC-order -1 -r GRCh37 -p /fdb/moc
ha/GRCh37/cnp.grch37.bed.gz -o output.vcf input.vcf.gz; Date=Mon Sep 28 09:51:10 2020
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  TCGA-B5-A0JV-01A-11D-A10B-09     TC
GA-B5-A0JV-10A-01W-A10C-09
1       3703493 .       G       A       .       .       .       GT:AD:DP:Ldev:Bdev:Bdev_Phase   0/1:57,9:.
:0:0:0  0/0:24,0:.:0:0:0
1       19433445        .       C       T       .       .       .       GT:AD:DP:Ldev:Bdev:Bdev_Phase    0/
1:63,13:.:0:0:0 0/0:115,0:.:0:0:0
1       21804793        .       C       T       .       .       .       GT:AD:DP:Ldev:Bdev:Bdev_Phase    0/
1:181,47:.:0:0:0        0/0:347,0:.:0:0:0
1       21952879        .       C       T       .       .       .       GT:AD:DP:Ldev:Bdev:Bdev_Phase    0/
...

```

Exit the mocha application:
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. mocha.sh). For example:



```

#!/bin/bash
set -e
module load mocha
export TMPDIR=/lscratch/$SLURM_JOBID
mocha  -g GRCh37 -p $MOCHA_REF/cnp.grch37.bed.gz -o output1.vcf input1.vcf.gz
mocha  -g GRCh37 -p $MOCHA_REF/cnp.grch37.bed.gz -o output2.vcf input2.vcf.gz
...

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=16 --mem=30g --gres=lscratch:20 mocha.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. mocha.swarm). For example:



```

mocha  -g GRCh37 -p $MOCHA_REF/cnp.grch37.bed.gz -o output1.vcf input1.vcf.gz
mocha  -g GRCh37 -p $MOCHA_REF/cnp.grch37.bed.gz -o output2.vcf input2.vcf.gz
...

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f mocha.swarm [-g 30] [-t 16] --module mocha
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module TEMPLATE Loads the TEMPLATE module for each subjob in the swarm 
 | |
 | |
 | |








