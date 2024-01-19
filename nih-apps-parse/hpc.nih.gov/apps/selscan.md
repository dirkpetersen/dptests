

document.querySelector('title').textContent = 'selscan: haplotype based scans for selection ';
**selscan: haplotype based scans for selection** 


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



selscan is a tool for haplotype-based scans to detect natural selection, 
which are useful to identify recent or ongoing positive selection in genomes.
It is an efficient multithreaded application that implements Extended Haplotype Homozygosity (EHH), 
Integrated Haplotype Score (iHS), and Cross-population EHH (XPEHH). 
selscan accepts phased genotypes in multiple formats, including TPED.



### References:


* Zachary A. Szpiech and Ryan D. Hernandez   

*selscan: An Efficient Multithreaded Program to Perform EHH-Based Scans for Positive Selection*  

 [Molecular Biology and Evolution](https://academic.oup.com/mbe/article/31/10/2824/1012603) 2014, **31**
(10):2824–2827 doi:10.1093/molbev/msu211


Documentation
* [selscan Github page](https://github.com/szpiech/selscan)
* [selscan User Manual](https://github.com/szpiech/selscan/blob/master/manual/selscan-manual.pdf)


Important Notes
* Module Name: selscan (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set
	+ **SELSCAN\_HOME**  installation directory
	+ **SELSCAN\_BIN**       executable directory
	+ **SELSCAN\_SRC**       source code directory
	+ **SELSCAN\_DATA**  sample data and checkpoints directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=10g --gres=lscratch:10 -c8**
[user@cn3101 ~]$**module load selscan** 
[+] Loading selscan  1.3.0    

```

The available executables are:

```

[user@cn3101]$ **ls $SELSCAN\_BIN** 
norm  selscan 

```

In particular, the command line options of the executable selscan are as follows: 

```

[user@cn3101]$ **selscan --help**
selscan v1.3.0

selscan v1.3.0 -- a program to calculate EHH-based scans for positive selection in genomes.
Source code and binaries can be found at .

selscan currently implements EHH, iHS, XP-EHH, and nSL.

Citations:

selscan: ZA Szpiech and RD Hernandez (2014) MBE 31: 2824-2827.
iHH12: R Torres et al. (2018) PLoS Genetics 15: e1007898.
 N Garud et al. (2015) PLoS Genetics 11: 1–32.
nSL: A Ferrer-Admetlla et al. (2014) MBE 31: 1275-1291.
XP-nSL: Szpiech et al. (2020) bioRxiv doi:
 https://doi.org/10.1101/2020.05.19.104380.
XP-EHH: PC Sabeti et al. (2007) Nature 449: 913–918.
 K Wagh et al. (2012) PloS ONE 7: e44751.
iHS: BF Voight et al. (2006) PLoS Biology 4: e72.
EHH: PC Sabeti et al. (2002) Nature 419: 832–837.

To calculate EHH:

./selscan --ehh  --vcf  --map  --out 

To calculate iHS:

./selscan --ihs --vcf  --map  --out 

To calculate nSL:

./selscan --nsl --vcf  --out 

To calculate XP-nSL:

./selscan --xpnsl --vcf  --vcf-ref  --out 

To calculate iHH12:

./selscan --ihh12 --vcf  --map  --out 

To calculate XP-EHH:

./selscan --xpehh --vcf  --vcf-ref  --map  --out 

----------Command Line Arguments----------

--alt : Set this flag to calculate homozygosity based on the sum of the
 squared haplotype frequencies in the observed data instead of using
 binomial coefficients.
 Default: false

--cutoff : The EHH decay cutoff.
 Default: 0.05

--ehh : Calculate EHH of the '1' and '0' haplotypes at the specified
 locus. Output:   <'1' EHH> <'0' EHH>
 Default: \_\_NO\_LOCUS\_\_

--ehh-win : When calculating EHH, this is the length of the window in bp
 in each direction from the query locus.
 Default: 100000

--gap-scale : Gap scale parameter in bp. If a gap is encountered between
 two snps > GAP\_SCALE and < MAX\_GAP, then the genetic distance is
 scaled by GAP\_SCALE/GAP.
 Default: 20000

--hap : A hapfile with one row per haplotype, and one column per
 variant. Variants should be coded 0/1
 Default: \_\_hapfile1

--help : Prints this help dialog.
 Default: false

--ihh12 : Set this flag to calculate iHH12.
 Default: false

--ihs : Set this flag to calculate iHS.
 Default: false

--ihs-detail : Set this flag to write out left and right iHH scores for '1' and '0' in addition to iHS.
 Default: false

--keep-low-freq : Include low frequency variants in the construction of your haplotypes.
 Default: false

--maf : If a site has a MAF below this value, the program will not use
 it as a core snp.
 Default: 0.05

--map : A mapfile with one row per variant site.
 Formatted    .
 Default: \_\_mapfile

--max-extend : The maximum distance an EHH decay curve is allowed to extend from the core.
 Set <= 0 for no restriction.
 Default: 1000000

--max-extend-nsl : The maximum distance an nSL haplotype is allowed to extend from the core.
 Set <= 0 for no restriction.
 Default: 100

--max-gap : Maximum allowed gap in bp between two snps.
 Default: 200000

--nsl : Set this flag to calculate nSL.
 Default: false

--out : The basename for all output files.
 Default: outfile

--pi : Set this flag to calculate mean pairwise sequence difference in a sliding window.
 Default: false

--pi-win : Sliding window size in bp for calculating pi.
 Default: 100

--pmap : Use physical map instead of a genetic map.
 Default: false

--ref : A hapfile with one row per haplotype, and one column per
 variant. Variants should be coded 0/1. This is the 'reference'
 population for XP-EHH calculations. Ignored otherwise.
 Default: \_\_hapfile2

--skip-low-freq : \*\*This flag is now on by default. If you want to include low frequency variants
in the construction of your haplotypes please use the --keep-low-freq flag.
 Default: false

--threads : The number of threads to spawn during the calculation.
 Partitions loci across threads.
 Default: 1

--tped : A TPED file containing haplotype and map data.
 Variants should be coded 0/1
 Default: \_\_hapfile1

--tped-ref : A TPED file containing haplotype and map data.
 Variants should be coded 0/1. This is the 'reference'
 population for XP-EHH calculations and should contain the same number
 of loci as the query population. Ignored otherwise.
 Default: \_\_hapfile2

--trunc-ok : If an EHH decay reaches the end of a sequence before reaching the cutoff,
 integrate the curve anyway (iHS and XPEHH only).
 Normal function is to disregard the score for that core.
 Default: false

--vcf : A VCF file containing haplotype data.
 A map file must be specified with --map.
 Default: \_\_hapfile1

--vcf-ref : A VCF file containing haplotype and map data.
 Variants should be coded 0/1. This is the 'reference'
 population for XP-EHH calculations and should contain the same number
 of loci as the query population. Ignored otherwise.
 Default: \_\_hapfile2

--wagh : Set this flag to calculate XP-EHH using definition of EHH which
 separates core SNP alleles in the denominator.
 Default: false

--xpehh : Set this flag to calculate XP-EHH.
 Default: false

--xpnsl : Set this flag to calculate XP-nSL.
 Default: false

```

To perform training of the predictor network using this executable, copy sample data to the current folder:

```

[user@cn3101]$ **cp $SELSCAN\_DATA/\* .**

```

A sample command to run selscan:

```

[user@cn3101]$ **selscan --ehh Locus1 --hap example2.pop2.hap --map example2.map --out my\_output**
selscan v1.3.0
Opening example2.pop2.hap...
Loading 125 haplotypes and 12920 loci...
Opening example2.map...
Loading map data for 12920 loci
Found Locus1 in data.
--skip-low-freq set. Removing all variants < 0.05.
Removed 8039 low frequency variants.

```

The command the following output files:

```

my_output.ehh.Locus1.log
my_output.ehh.Locus1.out
my_output.ehh.Locus1.out.anc.colormap
my_output.ehh.Locus1.out.der.colormap

```

End the interactive session:

```

[user@cn3101 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```





