

document.querySelector('title').textContent = 'coverageMaster: CNV detection and visualization from NGS short reads';
coverageMaster: CNV detection and visualization from NGS short reads


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



CoverageMaster (CoM) is a copy number variation (CNV) calling algorithm i
ased on depth-of-coverage maps designed to detect CNVs of any size i
n exome [whole exome sequencing (WES)] and genome [whole genome sequencing (WGS)] data. 
The core of the algorithm is the compression of sequencing coverage data 
in a multiscale Wavelet space and the analysis through an iterative Hidden Markov Model. 



### References:


* Melivoia Rapti, Yassine Zouaghi, Jenny Meylan, Emmanuelle Ranza, Stylianos E Antonarakis, Federico A Santoni   

*CoverageMaster: comprehensive CNV detection and visualization from NGS short reads for genetic medicine applications*    

[Briefings in Bioinformatics, Volume 23, Issue 2, March 2022, bbac049, https://doi.org/10.1093/bib/bbac049](https://academic.oup.com/bib/article/23/2/bbac049/6537346?login=true) Published: 26 February 2022


Documentation
* [coverageMaster github page](https://github.com/fredsanto/coverageMaster)


Important Notes
* Module Name: coverageMaster (see [the modules page](/apps/modules.html) for more information)
* Unusual environment variables set
	+ **COM\_HOME**  installation directory
	+ **COM\_BIN**  executable directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive -c 16 --mem 45g --gres=lscratch:20**
[user@cn3144 ~]$ **module load coverageMaster**
[+] Loading singularity  3.8.5-1  on cn3144
[+]  Loading coverageMaster 20220706  ...
[user@cn3144 ~]$  **coverageMaster.py -h** 
Usage: coverageMaster.py [options] <cov_file> <stats_file> <gene list(file or comma separated gene names)|region(chr:start-end)> -r <reference.cov> -o <output_px>

Options:
  -h, --help            show this help message and exit
  -c CONTROL, --control=CONTROL
                        <optional> txt file with a control file per line (.cov
                        with .report.txt in same folder)
  -s SINGLE_CONTROL, --single-control=SINGLE_CONTROL
                        <optional> single control file (.cov with .report.txt
                        in same folder)
  -r REF, --ref=REF     reference file
  -o OUTPUT_PX, --out=OUTPUT_PX
                        output prefix
  -g CGD, --cgd=CGD     <optional> clinical genomic database
  -f, --force           <optional> force output
  -e EXONS, --exons=EXONS
                         <optional> n. of extra exons
  -x OFFSET, --offset=OFFSET
                         <optional> offset to ref
  -d WID, --width=WID    <optional> d*std
  -l LEV, --level=LEV    <optional> max wavelet level
  -m MINLEV, --minlevel=MINLEV
                         <optional> min wavelet level
  -w, --wig              <optional> write wig
  -b, --bed              <optional> BED input
  -k DGV, --dgv=DGV      <optional> DGV file
[user@cn3144 ~]$ **git clone https://github.com/fredsanto/coverageMaster**
[user@cn3144 ~]$ **cd coverageMaster/DEMO**
[user@cn3144 ~]$ **gene=PGM1 && co=control.PGM1.cov && python ../coverageMaster.py test.PGM1.cov test.PGM1.report.txt $gene -s $co -r ref.PGM1 -o test.PGM1**
Fri Sep  9 04:04:48 2022 - --------------------------------------------------------------------------------
CoverageMaster is warming up
Fri Sep  9 04:04:48 2022 - Query index created
Fri Sep  9 04:04:48 2022 - Reference index created
Fri Sep  9 04:04:48 2022 - Processing Gene: XIST - Region chrX:73040485-73072588
chrX not in my list - skipped
/opt/conda/envs/coveragemaster/lib/python3.6/site-packages/numpy/core/fromnumeric.py:3118: RuntimeWarning: Mean of empty slice.
  out=out, **kwargs)
/opt/conda/envs/coveragemaster/lib/python3.6/site-packages/numpy/core/_methods.py:85: RuntimeWarning: invalid value encountered in double_scalars
  ret = ret.dtype.type(ret / rcount)
Fri Sep  9 04:04:48 2022 - Control index created
Fri Sep  9 04:04:48 2022 - Processing Gene: XIST - Region chrX:73040485-73072588
chrX not in my list - skipped
Fri Sep  9 04:04:48 2022 - Loop with control control.PGM1.cov
Fri Sep  9 04:04:48 2022 - Processing Gene: PGM1 - Region chr1:64058946-64125916
Fri Sep  9 04:04:48 2022 - Processing Gene: PGM1 - Region chr1:64058946-64125916
Fri Sep  9 04:04:48 2022 - PGM1:Zooming 4
Fri Sep  9 04:04:48 2022 - PGM1:Zooming 3
Fri Sep  9 04:04:48 2022 - PGM1:Zooming 2
Fri Sep  9 04:04:48 2022 - PGM1:Zooming 1
user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```





