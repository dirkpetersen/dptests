

document.querySelector('title').textContent = "Chromeister";
CHROMEISTER on Biowulf


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



Chromeister: An ultra fast, heuristic approach to detect conserved signals in extremely large pairwise genome comparisons.



### References:


* Pérez-Wohlfeil, E., Diaz-del-Pino, S. & Trelles, O.
 [**Ultra-fast genome comparison for large-scale genomic experiments.**](https://doi.org/10.1038/s41598-019-46773-w)
*Scientific reports 9, no. 1 (2019): 1-10.*


Documentation
* [Chromeister Main Site](https://github.com/estebanpw/chromeister)


Important Notes
* Module Name: chromeister (see [the modules page](/apps/modules.html) for more information)
 * Multithreaded
* Environment variables set 
	+ CHROMEISTER\_HOME
	+ SCRIPTS adds scripts from /usr/local/apps/chromeister/1.5.a/scripts to path



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive**
[user@cn4274 ~]$ **module load chromeister**
[+] Loading chromeister  1.5.a  on cn4274
[+] Loading singularity  3.10.5  on cn4274
[+] Loading gcc  11.3.0  ...
[+] Loading HDF5  1.12.2
[+] Loading netcdf  4.9.0
[-] Unloading gcc  11.3.0  ...
[+] Loading gcc  11.3.0  ...
[+] Loading openmpi/4.1.3/gcc-11.3.0  ...
[+] Loading pandoc  2.18  on cn4274
[+] Loading pcre2  10.40
[+] Loading R 4.3.0

```


Example
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Execute Chromeister binary



```

#copy over test-data
[user@cn4274 ~]$ **cp -a /usr/local/apps/chromeister/test-data .**
[user@cn4274 ~]$ **cd test-data**
#run chromeister with inputs
[user@cn4274 test-data]$ **CHROMEISTER -query mycoplasma-232.fasta \**
**-db mycoplasma-232.fasta \**
**-out mycoplasma-232-7422.mat \**
**-dimension 500 && Rscript ${SCRIPTS}/compute\_score.R mycoplasma-232-7422.mat 500**
[INFO] Generating a 500x500 matrix
[INFO] Loading database
99%...[INFO] Database loaded and of length 892758.
[INFO] Ratios: Q [1.785516e+03] D [1.785516e+03]. Lenghts: Q [892758] D [892758]
[INFO] Pixel size: Q [5.600622e-04] D [5.600622e-04].
[INFO] Computing absolute hit numbers.
99%...Scanning hits table.
99%...
[INFO] Query length 892758.
[INFO] Writing matrix.
[INFO] Found 25819 unique hits for z = 4.
0

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








