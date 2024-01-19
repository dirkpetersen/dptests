

document.querySelector('title').textContent = 'hichipper on Biowulf';
hichipper on Biowulf


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



The hichipper package implements our data processing and quality control pipeline for HiChIP data. This package takes output from a HiC-Pro run and a sample manifest file (.yaml) that coordinates optional high-quality peaks (identified through ChIP-Seq) and restriction fragment locations (see folder here) as input and produces output that can be used to 1) determine library quality, 2) identify and characterize DNA loops and 3) interactively visualize loops.



### References:


* Lareau, C.A. and Aryee, M.J. (2017) "hichipper: A preprocessing pipeline for assessing library quality and DNA loops from HiChIP data." bioRxiv doi: [10.1101/192302](https://doi.org/10.1101/192302)


Documentation
* [hichipper main site](http://hichipper.readthedocs.io/en/latest/)


Important Notes
* Module Name: hichipper (see [the modules page](/apps/modules.html) for more information)
* Example files in /fdb/hichipper/test



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --gres lscratch:1**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job
[user@cn3144 ~]$ **cd /lscratch/$SLURM\_JOB\_ID**
[user@cn3144 46116226]$ **module load hichipper**
[+] Loading samtools 1.5 ...
[+] Loading bedtools 2.25.0
[+] Loading pandoc 1.15.0.6 ...
[+] Loading GSL 2.2.1 ...
[+] Loading Graphviz v2.38.0 ...
[+] Loading gdal 2.0 ...
[+] Loading proj 4.9.2 ...
[+] Loading gcc 4.9.1 ...
[+] Loading openmpi 1.10.0 for GCC 4.9.1
[+] Loading tcl_tk 8.6.3
[-] Unloading pandoc 1.15.0.6 ...
[+] Loading pandoc 1.15.0.6 ...
[+] Loading Zlib 1.2.8 ...
[+] Loading Bzip2 1.0.6 ...
[+] Loading pcre 8.38 ...
[+] Loading liblzma 5.2.2 ...
[-] Unloading Zlib 1.2.8 ...
[+] Loading Zlib 1.2.8 ...
[-] Unloading liblzma 5.2.2 ...
[+] Loading liblzma 5.2.2 ...
[+] Loading libjpeg-turbo 1.5.1 ...
[+] Loading tiff 4.0.7 ...
[+] Loading curl 7.46.0 ...
[+] Loading boost libraries v1.65 ...
[+] Loading R 3.4.0 on cn3108
[+] Loading hichipper, version 0.7.0...
[user@cn3144 46116226]$ **cp -r /fdb/hichipper/tests/ .**
[user@cn3144 46116226]$ **cd tests**
[user@cn3144 tests]$ **hichipper --out output1 yaml/one.yaml**
Thu Oct 19 12:17:16 EDT 2017: Starting hichipper pipeline v0.7.0
Thu Oct 19 12:17:16 EDT 2017: Parsing user parameters
Anchors removed due to excessive size (likely at ends of chromosomes): 0 
Thu Oct 19 12:17:20 EDT 2017: Processing d1
Thu Oct 19 12:17:20 EDT 2017: Total_PETs=164332898
Thu Oct 19 12:17:20 EDT 2017: Mapped_unique_quality_pairs=1888170
Thu Oct 19 12:17:20 EDT 2017: Mapped_unique_quality_valid_pairs=627403
Thu Oct 19 12:17:20 EDT 2017: Intersecting PETs with anchors
Thu Oct 19 12:17:20 EDT 2017: Finished the anchor merging.
Thu Oct 19 12:17:52 EDT 2017: Intrachromosomal_valid_small=163462
Thu Oct 19 12:17:53 EDT 2017: Intrachromosomal_valid_med=394048
Thu Oct 19 12:17:54 EDT 2017: Intrachromosomal_valid_large=69902
Thu Oct 19 12:17:54 EDT 2017: Total number of anchors used: 926
Thu Oct 19 12:17:54 EDT 2017: Total number of reads in anchors: 410414
Thu Oct 19 12:17:59 EDT 2017: Mapped_unique_intra_quality_anchor=20288
Thu Oct 19 12:17:59 EDT 2017: Mapped_unique_intra_quality_anchor_small=13632
Thu Oct 19 12:17:59 EDT 2017: Mapped_unique_intra_quality_anchor_med=6090
Thu Oct 19 12:17:59 EDT 2017: Mapped_unique_intra_quality_anchor_large=566
Thu Oct 19 12:17:59 EDT 2017: Loop_PETs=6090
Thu Oct 19 12:17:59 EDT 2017: 
['Rscript', u'output1/qcReport.R', '/usr/local/Anaconda/envs_app/hichipper/0.7.0/lib/python2.7/site-packages/hichipper', 'output1', '/lscratch/51844726/tests', '0.7.0', 'd1']
/lscratch/51844726/tests

processing file: qcReport_make.Rmd
  |...                                                              |   5%
   inline R code fragments

  |.......                                                          |  11%
label: unnamed-chunk-1 (with options) 
List of 3
 $ echo   : logi FALSE
 $ message: logi TRUE
 $ warning: logi TRUE

  |..........                                                       |  16%
   inline R code fragments

  |..............                                                   |  21%
label: unnamed-chunk-2 (with options) 
List of 3
 $ echo   : logi FALSE
 $ message: logi FALSE
 $ warning: logi FALSE

  |.................                                                |  26%
  ordinary text without R code

  |.....................                                            |  32%
label: unnamed-chunk-3 (with options) 
List of 4
 $ echo   : logi FALSE
 $ message: logi FALSE
 $ warning: logi FALSE
 $ results: chr "asis"

  |........................                                         |  37%
  ordinary text without R code

  |...........................                                      |  42%
label: unnamed-chunk-4 (with options) 
List of 4
 $ echo     : logi FALSE
 $ message  : logi FALSE
 $ warning  : logi FALSE
 $ out.width: chr "\\textwidth"

  |...............................                                  |  47%
  ordinary text without R code

  |..................................                               |  53%
label: unnamed-chunk-5 (with options) 
List of 4
 $ echo     : logi FALSE
 $ message  : logi FALSE
 $ warning  : logi FALSE
 $ out.width: chr "\\textwidth"

  |......................................                           |  58%
  ordinary text without R code

  |.........................................                        |  63%
label: unnamed-chunk-6 (with options) 
List of 6
 $ echo     : logi FALSE
 $ message  : logi FALSE
 $ warning  : logi FALSE
 $ results  : chr "asis"
 $ out.width: chr "\\textwidth"
 $ fig.width: num 7

  |............................................                     |  68%
  ordinary text without R code

  |................................................                 |  74%
label: unnamed-chunk-7 (with options) 
List of 4
 $ echo     : logi FALSE
 $ message  : logi FALSE
 $ warning  : logi FALSE
 $ out.width: chr "\\textwidth"

Using count as value column: use value.var to override.
  |...................................................              |  79%
  ordinary text without R code

  |.......................................................          |  84%
label: unnamed-chunk-8 (with options) 
List of 4
 $ echo     : logi FALSE
 $ message  : logi FALSE
 $ warning  : logi FALSE
 $ out.width: chr "\\textwidth"

  |..........................................................       |  89%
  ordinary text without R code

  |..............................................................   |  95%
label: unnamed-chunk-9 (with options) 
List of 4
 $ echo     : logi FALSE
 $ message  : logi FALSE
 $ warning  : logi FALSE
 $ out.width: chr "\\textwidth"

  |.................................................................| 100%
  ordinary text without R code


output file: qcReport_make.knit.md

/usr/local/apps/pandoc/1.15.0.6/bin/pandoc +RTS -K512m -RTS qcReport_make.utf8.md --to html --from markdown+autolink_bare_uris+ascii_identifiers+tex_math_single_backslash --output output1.hichipper.qcreport.html --smart --email-obfuscation none --self-contained --standalone --section-divs --template /usr/local/apps/R/gcc_4.9.1/library/rmarkdown/rmd/h/default.html --highlight-style pygments --css /usr/local/apps/R/gcc_4.9.1/library/rmarkdown/rmarkdown/templates/html_vignette/resources/vignette.css --include-in-header /lscratch/51844726/RtmpopswEW/rmarkdown-str60995cfb604b.html --mathjax --variable 'mathjax-url:https://mathjax.rstudio.com/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML' 

Output created: output1.hichipper.qcreport.html
Warning messages:
1: Transformation introduced infinite values in continuous x-axis 
2: Removed 13419 rows containing non-finite values (stat_bin). 
Processing: d1
[user@cn3144 tests]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. hichipper.sh). For example:



```

#!/bin/sh
set -e
module load hichipper

cp -r /fdb/hichipper/tests/ .
cd tests
hichipper --out output1 yaml/one.yaml

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] hichipper.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. hichipper.swarm). For example:



```

hichipper --out sample1 sample1.yaml
hichipper --out sample2 sample2.yaml
hichipper --out sample3 sample3.yaml
hichipper --out sample4 sample4.yaml

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f hichipper.swarm [-g #] [-t #] --module hichipper
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g #  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t #  Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module hichipper  Loads the hichipper module for each subjob in the swarm 
 | |
 | |
 | |








