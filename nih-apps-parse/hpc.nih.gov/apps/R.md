

document.querySelector('title').textContent = 'R/Bioconductor on Biowulf';




 .hl { background-color: #ffff99; }
 code {padding: 1px; background-color: #eeeeee; font-family: monospace;}
 dt {font-weight: bold; margin-top: 5px;}
 dd {padding-left: 2px; border-left: 1px solid #bbbbbb;}
 .btt {border: 1px solid silver;
 background-color: white;
 padding: 5px;
 position: relative;
 margin: 5px 0px 10px 10px;
 float: right;
 top: -25px;
 left: 10px;
 }

R/Bioconductor on Biowulf


|  |
| --- |
| 
Quick Links
[Changelog](#changes)
[Common pitfalls](#gotcha)
[Interactive R](#int)
[Installed packages](#packages)
[Manage your own packages](#install)
[R batch job](#sbatch)
[Swarm of R jobs](#swarm) 
[Rswarm](#rswarm)
[Using the parallel package](#parallel)
[Using the BiocParallel package](#biocparallel)
[Implicit multithreading](#threading)
[R MPI jobs](#rmpi)
[Run a shiny app on biowulf](#shiny)
[Notes for individual packages](#pkgnotes)
[Documentation](#doc)
 |



[![R logo](/images/Rlogo.jpg)](http://cran.r-project.org)
R is a language and environment for statistical computing and graphics. It
can be considered an open source decendant of the S language which was
developed by Chambers and colleagues at Bell Laboratories in the 1970s.


R is highly extensible and provides a wide variety of modern statistical
analysis methods combined with excellent graphical visualization capabilities
embedded in a programming language that supports procedural, functuional, and
object oriented programming styles. R natively provides operators for calculations
on arrays and matrices.


While individual single threaded R code is not expected to run any faster in
an [interactive session](/docs/userguide.html#int) on a compute node
than it would run on a modern desktop, Biowulf allows users to run many R jobs
concurrently or take advantage of the speedup provided by parallelizing R code
to a much greater degree than possible on a single workstation.


On biowulf, R modules are available for the minor releases (e.g. 4.2) which
will contain the newest patch level releases (e.g. 4.2.3).


 

Changelog
[top](#top)

June 23 2023: default location for `$R_LIBS_USER`
changed to `/data/$USER/R/rhel8/%v` with migration to RHEL8
May 2023: R/4.3.0 becomes the default R installation.
User visible changes:
 * Calling `&&` or `||` with LHS or (if evaluated) 
 RHS of length greater than one is now always an error, with a report of the form
 `'length = 4' in coercion to 'logical(1)'`. R 4.2.0 introduced
 a warning for this usage of conditional operators.
* Environment variable \_R\_CHECK\_LENGTH\_1\_LOGIC2\_ no longer has any effect.


 See [R NEWS](https://cran.r-project.org/doc/manuals/r-release/NEWS.html) for
 full details.
 
Nov 2022: R/4.2.2 becomes the default R installation.
Jul 2022: R/4.2.0 becomes the default R installation.
Some notable changes:
 * Calling `&&` or `||` with either argument of length greater
 than one now gives a warning (which it is intended will become
 an error).
* Calling if() or while() with a condition of length greater
 than one gives an error rather than a warning. Consequently,
 environment variable \_R\_CHECK\_LENGTH\_1\_CONDITION\_ no longer has
 any effect.
* The graphics engine version, R\_GE\_version, has been bumped to 15.



Jul 2021: R/4.1.0 becomes the default R installation.
Apr 2021: R/4.0.5 becomes the default R installation. OpenMPI is now version 4
Nov 2020: R/4.0.3 becomes the default R installation
Jun 2020: R/4.0.0 becomes the default R installation
For details see [R NEWS](http://cran.r-project.org/doc/manuals/r-release/NEWS.html). Notable changes:
 * R now uses `stringsAsFactors = FALSE` as the default
* There is a new syntax for specifying raw character constants similar to the one used in C++: r"(...)"


 As usual, many packages are pre-installed and private packages
 need to be re-installed.
 
 
Jun 2020: R/4.0.0 becomes the default R installation
Apr 2020: R/3.6.3 becomes the default R installation
Dec 2019: R/3.6.1 becomes the default R installation and R is now compiled with gcc 9.2.0
Jun 2019: R/3.6.0 becomes the default R installation
Feb 2019: R/3.5.2 becomes the default R installation
Jun 2018: Cluster update from RHEL6 to RHEL7

* R sessions are not allowed on helix any more.
* R is compiled against MKL for improved performance





Common pitfalls
[top](#top)

Implicit multithreading

 R can make use of [implicit multithreading](#threading) via two
 different mechanisms. One of them is regulated by the `OMP_NUM_THREADS`
 environment variable which is set to 1 by the R modules because leaving this
 variable unset can lead to R using as many threads as there are CPUs on a compute
 node thus overloading jobs. If you know your code can make effective use of those
 threads you can explicitly set `OMP_NUM_THREADS` to greated than 1 after
 loading the module. However, only a subset of code will be able to take advantage of
 this - don't expect an automatic speed increase.
 
`parallel::detectCores()` always detects all CPUs on a node

 R using one of the parallel packages (parallel, doParallel, ...) often overload
 their job allocation because they are using the `detectCores()` function
 from the parallel package to determine how many worker processes to use. However,
 this function returns the number of physical CPUs on a compute node irrespective of
 how many have been allocated to a job. Therefore, if not all CPUs are allocated to a job
 the job will be overloaded and perform poorly. See the section on the 
 [parallel](#parallel) package for more detail.
 
`BiocParallel` by default tries to use most CPUs on a node

`BiocParallel` is not aware of slurm and by default tries to use most
 of the CPUs on a node irrespetive of the slurm allocation. This can lead to overloaded
 jobs. See the section on the
 [BiocParallel](#biocparallel) package for more information on how to avoid
 this.
 
Poor scaling of parallel code

 Don't assume that you should allocate as many CPUs as possible to a parallel workload.
 Parallel efficiency often drops and in some cases allocating more CPUs may actually
 extend runtimes. If you use/implement parallel algorithms please measure scaling
 before submitting large numbers of such jobs.
 
Can't install packages in my private library

 R will attempt to install packages to `/data/R/rhel8/%v` (RHEL8) where `%v` is the two
 digit version of R (e.g. 4.1) or in the path set by `$R_LIBS_USER`. 
 However, R won't always automatically create that directory and in
 its absence will try to install to the central packge library which will fail. If you encounter
 installation failures please make sure the library directory for your version of R exists.
 
AnnotationHub or ExperimentHub error `No internet connection`

 The AnnotationHub and ExperimentHub packages and packages depending on them need to connect to the
 internet via a proxy. When using AnnotationHub or ExperimentHub directly, a proxy can be specified
 explicitly in the call to set up the Hub. However, if they are used indirectly that is not possible.
 Instead, define the proxy either using environment variables EXPERIMENT\_HUB\_PROXY/ANNOTATION\_HUB\_PROXY
 or by setting options in R with `setAnnotationHubOption("PROXY", Sys.getenv("http_proxy"))`
 or the corresponding `setExperimentHubOption` function.
 
Updating broken packages installed in home directory

 When R and/or the centrally installed R packages are updated, packages installed in your
 private library may break or result in other packages not loading. The most common error
 results from a locally installed `rlang` package. Look for errors that
 include your private library path in the errors. All locally
 installed package can be updated with
 
```

> my.lib <- .libPaths()[1]  # check first that .libPaths()[1] is indeed the path to your library
> my.pkgs <- list.files(my.lib)
> library(pacman)
> p_install(pkgs, character.only=T, lib=my.lib)

```

 An easy way to fix this error is also to delete the locally installed rlang with
 
```

$ rm -rf ~/R/4.2/library/rlang  # replace 4.2 with the R major.minor version you are using
        
```

Re-install packages to 1) data directory from home directory and/or 2) different R version

 Since R\_LIBS\_USER was relocated to data directory, the R packages installed in your
 private library (at ~/R/) need to be reinstalled. And/or if you need to update to a newer R version by reinstalling
 all the packages from an older version of R. Let's create the private library directory for the new version of R first (e.g. R/4.3):

```

$ mkdir -p /data/$USER/R/rhel8/4.3/
        
```

For example, for packages installed under R/4.2, you can re-install them by creating a list of installed libraries, find 
the ones that are not yet installed under data directory or not the same R version you are using, then re-install them with (please replace apptest with your username):

```

> packages<-installed.packages(lib.loc="/home/apptest/R/4.2/library")[,"Package"]
> toInstall<-setdiff(packages,installed.packages(loc.lib="/data/apptest/R/rhel8/4.3/")[,"Package"])
> BiocManager::install(toInstall)

```



R will automatically use [lscratch](https://hpc.nih.gov/docs/userguide.html#local)
for temporary files if it has been allocated. Therefore we highly recommend users always
allocate a minimal amount of lscratch of 1GB plus whatever lscratch storage is required by
your code.



Interactive R
[top](#top)
Allocate an [interactive session](/docs/userguide.html#int) for
interactive R work. Note that R sessions are not allowed on the login node
nor helix.



```

[user@biowulf]$ **sinteractive --gres=lscratch:5**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

```

There may be multiple versions of R available. An easy way of selecting the
version is to use [modules](/apps/modules.html).To see the modules
available, type



```

[user@cn3144 ~]$ **module -r avail '^R$'**

--------------- /usr/local/lmod/modulefiles ------------------
R/3.4    R/3.4.3    R/3.4.4    R/3.5 (D)    R/3.5.0    R/3.5.2

```

Set up your environment and start up an R session



```

[user@cn3144 ~]$ **module load R/3.5**
[user@cn3144 ~]$ **R**
R version 3.5.0 (2018-04-23) -- "Joy in Playing"
Copyright (C) 2018 The R Foundation for Statistical Computing
Platform: x86_64-pc-linux-gnu (64-bit)

R is free software and comes with ABSOLUTELY NO WARRANTY.
You are welcome to redistribute it under certain conditions.
Type 'license()' or 'licence()' for distribution details.

  Natural language support but running in an English locale

R is a collaborative project with many contributors.
Type 'contributors()' for more information and
'citation()' on how to cite R or R packages in publications.

Type 'demo()' for some demos, 'help()' for on-line help, or
'help.start()' for an HTML browser interface to help.
Type 'q()' to quit R.

> library(tidyverse)
── Attaching packages ─────────────────────────────────────── tidyverse 1.2.1 ──
✔ ggplot2 2.2.1     ✔ purrr   0.2.4
✔ tibble  1.4.2     ✔ dplyr   0.7.4
✔ tidyr   0.8.0     ✔ stringr 1.3.0
✔ readr   1.1.1     ✔ forcats 0.3.0
── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
✖ dplyr::filter() masks stats::filter()
✖ dplyr::lag()    masks stats::lag()
> [...lots of work...]
> q()
Save workspace image? [y/n/c]: n

```

A rudimentary graphical interface is available if the sinteractive
session was started from a session with X11 forwarding enabled:



```

[user@cn3144 ~]$ **R --gui=Tk**

```

However, [RStudio](RStudio.html) is a much better interface
with many advanced features.


Don't forget to exit the interactive session



```

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Installed packages
[top](#top)
Packages installed in the current default R environment






| Package | Version |
| --- | --- |
| spatial | 7.3-16 |
| viridisLite | 0.4.2 |
| urltools | 1.7.3 |
| trapezoid | 2.0-2 |
| plotly | 4.10.2 |
| oligoClasses | 1.62.0 |
| numDeriv | 2016.8-1.1 |
| methylumi | 2.46.0 |
| mcmc | 0.9-7 |
| matrixcalc | 1.0-6 |
| latticeExtra | 0.6-30 |
| kdecopula | 0.9.2 |
| ipred | 0.9-14 |
| GWASExactHW | 1.01 |
| groHMM | 1.34.0 |
| GOSemSim | 2.26.1 |
| geosphere | 1.5-18 |
| DESeq2 | 1.40.2 |
| DEoptimR | 1.1-3 |
| cn.mops | 1.46.0 |
| brio | 1.1.3 |
| AnnotationForge | 1.42.2 |
| TxDb.Celegans.UCSC.ce6.ensGene | 3.2.2 |
| strucchange | 1.5-3 |
| stabledist | 0.7-1 |
| shadowtext | 0.1.2 |
| rmarkdown | 2.25 |
| RcppTOML | 0.2.2 |
| ProtGenerics | 1.32.0 |
| PMCMR | 4.4 |
| MultiAssayExperiment | 1.26.0 |
| msm | 1.7 |
| MLEcens | 0.1-7 |
| magic | 1.6-1 |
| JavaGD | 0.6-5 |
| ENmix | 1.36.08 |
| distillery | 1.2-1 |
| cpp11 | 0.4.6 |
| conflicted | 1.2.0 |
| caTools | 1.18.2 |
| bsplus | 0.1.4 |
| blob | 1.2.4 |
| ALL | 1.42.0 |
| ADGofTest | 0.3 |
| statmod | 1.5.0 |
| rcmdcheck | 1.4.0 |
| ks | 1.14.1 |
| keras | 2.13.0 |
| httr | 1.4.7 |
| Heatplus | 3.8.0 |
| europepmc | 0.4.3 |
| data.table | 1.14.8 |
| Cubist | 0.4.2.1 |
| CNEr | 1.36.0 |
| clusterRepro | 0.9 |
| bibtex | 0.5.1 |
| xfun | 0.40 |
| unisensR | 0.3.3 |
| spam | 2.9-1 |
| slingshot | 2.8.0 |
| NOISeq | 2.44.0 |
| Lmoments | 1.3-1 |
| interval | 1.1-1.0 |
| gfonts | 0.2.0 |
| DBI | 1.1.3 |
| ChIPpeakAnno | 3.34.1 |
| bindata | 0.9-20 |
| beadarray | 2.50.0 |
| RRPP | 1.4.0 |
| Rfast | 2.0.8 |
| PhysicalActivity | 0.2-4 |
| pheatmap | 1.0.12 |
| MplusAutomation | 1.1.0 |
| lifecycle | 1.0.3 |
| lars | 1.3 |
| Haplin | 7.3.0 |
| DaMiRseq | 2.12.0 |
| commonmark | 1.9.0 |
| ade4 | 1.7-22 |
| parallel | 4.3.0 |
| tfdatasets | 2.9.0 |
| remotes | 2.4.2.1 |
| plyr | 1.8.9 |
| org.Hs.eg.db | 3.17.0 |
| mize | 0.2.4 |
| kpeaks | 1.1.0 |
| flux | 0.3-0.1 |
| DAAG | 1.25.4 |
| bindr | 0.1.1 |
| Bhat | 0.9-12 |
| vipor | 0.4.5 |
| tximeta | 1.18.3 |
| tkWidgets | 1.78.0 |
| sciplot | 1.2-0 |
| refund | 0.1-32 |
| reactome.db | 1.84.0 |
| libcoin | 1.0-10 |
| cn.farms | 1.48.0 |
| cli | 3.6.1 |
| BiocSingular | 1.16.0 |
| bayestestR | 0.13.1 |
| slider | 0.3.1 |
| shinydashboard | 0.7.2 |
| sampling | 2.9 |
| RMTstat | 0.3.1 |
| raster | 3.6-26 |
| pls | 2.8-2 |
| miniUI | 0.1.1.1 |
| gsubfn | 0.7 |
| grr | 0.9.5 |
| GOstats | 2.66.0 |
| GENIE3 | 1.22.0 |
| gdata | 3.0.0 |
| future.apply | 1.11.0 |
| edgeR | 3.42.4 |
| dfoptim | 2023.1.0 |
| ChIPseeker | 1.36.0 |
| ccaPP | 0.3.3 |
| class | 7.3-22 |
| tensorA | 0.36.2 |
| spdep | 1.2-8 |
| RNifti | 1.5.0 |
| pROC | 1.18.4 |
| PBSmodelling | 2.68.8 |
| partykit | 1.2-20 |
| mstate | 0.3.2 |
| knn.covertree | 1.0 |
| GPArotation | 2023.8-1 |
| ggformula | 0.10.4 |
| extraDistr | 1.9.1 |
| DMRcatedata | 2.18.0 |
| cmprsk | 2.2-11 |
| batchelor | 1.16.0 |
| batch | 1.1-5 |
| additivityTests | 1.1-4.1 |
| svd | 0.5.4.1 |
| snow | 0.4-4 |
| Publish | 2023.01.17 |
| plotmo | 3.6.2 |
| plgem | 1.72.0 |
| parsedate | 1.3.1 |
| parathyroidSE | 1.38.0 |
| manhattanly | 0.3.0 |
| leafcutter | 0.2.9 |
| DRIMSeq | 1.28.0 |
| diffusionMap | 1.2.0 |
| bit | 4.0.5 |
| zlibbioc | 1.46.0 |
| rvinecopulib | 0.6.3.1.1 |
| qvcalc | 1.0.3 |
| pingr | 2.0.2 |
| pbdMPI | 0.4-6 |
| oro.nifti | 0.11.4 |
| mixsqp | 0.3-48 |
| MAGeCKFlute | 2.4.0 |
| e1071 | 1.7-13 |
| ashr | 2.2-63 |
| arsenal | 3.6.3 |
| TeachingDemos | 2.12 |
| R.matlab | 3.7.0 |
| rle | 0.9.2 |
| R2WinBUGS | 2.1-21 |
| markdown | 1.11 |
| ggm | 2.5 |
| fansi | 1.0.5 |
| clue | 0.3-65 |
| CHsharp | 0.4 |
| timeSeries | 4031.107 |
| sva | 3.48.0 |
| SummarizedExperiment | 1.30.2 |
| pec | 2023.04.12 |
| MatchIt | 4.5.5 |
| ICsurv | 1.0.1 |
| graph | 1.78.0 |
| ggthemes | 4.2.4 |
| ggdendro | 0.1.23 |
| ggalluvial | 0.12.5 |
| genetics | 1.3.8.1.3 |
| fontBitstreamVera | 0.1.1 |
| DNAcopy | 1.74.1 |
| discretization | 1.0-1.1 |
| colorRamps | 2.3.1 |
| ATACseqQC | 1.24.0 |
| tiff | 0.1-11 |
| survJamda.data | 1.0.2 |
| Rmpfr | 0.9-3 |
| rlecuyer | 0.3-7 |
| rafalib | 1.0.0 |
| pd.hugene.2.0.st | 3.14.1 |
| mygene | 1.36.0 |
| msa | 1.32.0 |
| geonames | 0.999 |
| FSelector | 0.34 |
| fracdiff | 1.5-2 |
| ELMER.data | 2.24.0 |
| EGSEA | 1.28.0 |
| dnet | 1.1.7 |
| CODEX | 1.32.0 |
| basilisk.utils | 1.12.1 |
| KernSmooth | 2.23-21 |
| wateRmelon | 2.6.0 |
| Surrogate | 3.2.1 |
| ROTS | 1.28.0 |
| questionr | 0.7.8 |
| pkgload | 1.3.3 |
| pkgbuild | 1.4.2 |
| optigrab | 0.9.2.1 |
| logcondens | 2.1.8 |
| JADE | 2.0-4 |
| inaparc | 1.2.0 |
| ids | 1.0.1 |
| gsl | 2.1-8 |
| ggtree | 3.8.2 |
| fdrtool | 1.2.17 |
| blockmodeling | 1.1.5 |
| utils | 4.3.0 |
| graphics | 4.3.0 |
| zCompositions | 1.4.1 |
| vioplot | 0.4.0 |
| rbenchmark | 1.0.0 |
| moments | 0.14.1 |
| locfit | 1.5-9.8 |
| inline | 0.3.19 |
| idr | 1.3 |
| hgu133plus2probe | 2.18.0 |
| gld | 2.6.6 |
| GGIRread | 0.3.1 |
| enrichR | 3.2 |
| dropbead | 0.3.1 |
| CCA | 1.2.2 |
| MASS | 7.3-60 |
| grDevices | 4.3.0 |
| xopen | 1.0.0 |
| tidyselect | 1.2.0 |
| seqLogo | 1.66.0 |
| rvest | 1.0.3 |
| robustbase | 0.99-0 |
| rjson | 0.2.21 |
| qs | 0.25.5 |
| praise | 1.0.0 |
| KMsurv | 0.1-5 |
| INPower | 1.36.0 |
| HilbertVis | 1.58.0 |
| hexView | 0.3-4 |
| geometry | 0.4.7 |
| EBSeqHMM | 1.34.0 |
| DMRcate | 2.14.1 |
| corpcor | 1.6.10 |
| spData | 2.3.0 |
| rotl | 3.1.0 |
| rex | 1.2.1 |
| RCircos | 1.2.2 |
| rbokeh | 0.5.2 |
| qpdf | 1.3.2 |
| psych | 2.3.9 |
| progress | 1.2.2 |
| PBSmapping | 2.73.2 |
| mr.raps | 0.2 |
| miscTools | 0.6-28 |
| MassSpecWavelet | 1.66.0 |
| LSD | 4.1-0 |
| gss | 2.2-7 |
| batchtools | 0.9.17 |
| airway | 1.20.0 |
| survival | 3.5-5 |
| TH.data | 1.1-2 |
| SPAtest | 3.1.2 |
| randtoolbox | 2.0.4 |
| ncbit | 2013.03.29.1 |
| hexbin | 1.28.3 |
| dagitty | 0.3-1 |
| bsseq | 1.36.0 |
| bslib | 0.5.1 |
| apeglm | 1.22.1 |
| spatstat.random | 3.1-6 |
| ReactomePA | 1.44.0 |
| pctGCdata | 0.3.0 |
| outliers | 0.15 |
| mlmRev | 1.0-8 |
| inum | 1.0-5 |
| illuminaio | 0.42.0 |
| hopach | 2.60.0 |
| exomeCopy | 1.46.0 |
| densityClust | 0.3.2 |
| bit64 | 4.0.5 |
| annaffy | 1.72.0 |
| tidygraph | 1.2.3 |
| rpf | 1.0.14 |
| rGADEM | 2.48.0 |
| plier | 1.70.0 |
| httr2 | 0.2.3 |
| git2r | 0.32.0 |
| gargle | 1.5.2 |
| gage | 2.50.0 |
| EnhancedVolcano | 1.18.0 |
| dqrng | 0.3.1 |
| ada | 2.0-5 |
| threejs | 0.3.3 |
| spacetime | 1.3-0 |
| Ringo | 1.64.0 |
| ncdf4 | 1.21 |
| MEGENA | 1.3.7 |
| ensemblVEP | 1.42.0 |
| DropletUtils | 1.20.0 |
| aroma.core | 3.3.0 |
| snpStats | 1.50.0 |
| semTools | 0.5-6 |
| rticles | 0.25 |
| riskRegression | 2023.09.08 |
| ResidualMatrix | 1.10.0 |
| R.cache | 0.16.0 |
| pathview | 1.40.0 |
| optmatch | 0.10.6 |
| motifmatchr | 1.22.0 |
| maSigPro | 1.72.0 |
| ipw | 1.2 |
| DDRTree | 0.1.5 |
| cqn | 1.46.0 |
| contrast | 0.24.2 |
| AUCell | 1.22.0 |
| aggregation | 1.0.1 |
| VGAM | 1.1-9 |
| vcd | 1.4-11 |
| styler | 1.10.2 |
| R.rsp | 0.45.0 |
| Rook | 1.2 |
| pwr | 1.3-0 |
| prediction | 0.3.14 |
| MotifDb | 1.42.0 |
| MAST | 1.26.0 |
| mapproj | 1.2.11 |
| lsmeans | 2.30-0 |
| Kendall | 2.2.1 |
| HDInterval | 0.2.4 |
| googlesheets4 | 1.1.1 |
| golubEsets | 1.42.0 |
| glmnet | 4.1-8 |
| futile.logger | 1.4.3 |
| boot | 1.3-28.1 |
| webdriver | 1.0.6 |
| tidybayes | 3.0.6 |
| RCurl | 1.98-1.12 |
| randtests | 1.0.1 |
| mipfp | 3.2.1 |
| HTqPCR | 1.54.0 |
| haven | 2.5.3 |
| ComplexHeatmap | 2.16.0 |
| tinytex | 0.48 |
| spatstat.utils | 3.0-3 |
| plm | 2.6-3 |
| pd.hg.u133.plus.2 | 3.12.0 |
| pbkrtest | 0.5.2 |
| PADOG | 1.42.0 |
| iterators | 1.0.14 |
| ica | 1.0-3 |
| hugene20sttranscriptcluster.db | 8.8.0 |
| gson | 0.1.0 |
| ggfortify | 0.4.16 |
| GetoptLong | 1.0.5 |
| dichromat | 2.0-0.1 |
| warp | 0.2.0 |
| rJava | 1.0-6 |
| pbdZMQ | 0.3-10 |
| pasilla | 1.28.0 |
| openssl | 2.1.1 |
| HTMLUtils | 0.1.8 |
| Gviz | 1.44.2 |
| ggbeeswarm | 0.7.2 |
| deconstructSigs | 1.8.0 |
| csaw | 1.34.0 |
| config | 0.3.2 |
| colourpicker | 1.3.0 |
| brglm | 0.7.2 |
| ballgown | 2.32.0 |
| datasets | 4.3.0 |
| remaCor | 0.0.16 |
| isva | 1.9 |
| isotone | 1.1-1 |
| filelock | 1.0.2 |
| ChIPsim | 1.54.0 |
| bio3d | 2.4-4 |
| base64enc | 0.1-3 |
| awsMethods | 1.1-1 |
| akima | 0.6-3.4 |
| signal | 0.7-7 |
| RWekajars | 3.9.3-2 |
| rngtools | 1.5.2 |
| RMySQL | 0.10.26 |
| Rcgmin | 2022-4.30 |
| mvQuad | 1.0-8 |
| mouse4302.db | 3.13.0 |
| marray | 1.78.0 |
| lobstr | 1.1.2 |
| fontawesome | 0.5.2 |
| deldir | 1.0-9 |
| ape | 5.7-1 |
| whisker | 0.4.1 |
| variancePartition | 1.30.2 |
| tweenr | 2.0.2 |
| SeqVarTools | 1.38.0 |
| random | 0.2.6 |
| profvis | 0.3.8 |
| prettyunits | 1.2.0 |
| PMCMRplus | 1.9.8 |
| openxlsx | 4.2.5.2 |
| MatrixEQTL | 2.3 |
| hunspell | 3.0.3 |
| hthgu133acdf | 2.18.0 |
| gridtext | 0.1.5 |
| ggraph | 2.1.0 |
| gamlss.data | 6.0-2 |
| flashClust | 1.01-2 |
| downloader | 0.4 |
| deepSNV | 1.46.0 |
| TitanCNA | 1.38.0 |
| sqldf | 0.4-11 |
| seriation | 1.5.1 |
| SCATE | 1.10.0 |
| rstan | 2.32.3 |
| loo | 2.6.0 |
| ggdist | 3.3.0 |
| genbankr | 1.27.0 |
| dr | 3.0.10 |
| DiceKriging | 1.6.0 |
| testthat | 3.2.0 |
| safe | 3.40.1 |
| rslurm | 0.6.2 |
| readstata13 | 0.10.1 |
| pryr | 0.1.6 |
| phyloseq | 1.44.0 |
| optimx | 2023-8.13 |
| modelr | 0.1.11 |
| mathjaxr | 1.6-0 |
| hu6800probe | 2.18.0 |
| diffHic | 1.32.0 |
| cosinor2 | 0.2.1 |
| coin | 1.4-3 |
| binom | 1.1-1.1 |
| tzdb | 0.4.0 |
| RUnit | 0.4.32 |
| rrcov | 1.7-4 |
| RPMM | 1.25 |
| rgdal | 1.6-6 |
| rematch2 | 2.1.2 |
| R6 | 2.5.1 |
| qvalue | 2.32.0 |
| multidplyr | 0.1.3 |
| laeken | 0.5.2 |
| invgamma | 1.1 |
| ijtiff | 2.3.3 |
| hdrcde | 3.4 |
| expint | 0.1-8 |
| dir.expiry | 1.8.0 |
| DiffBind | 3.10.1 |
| DEXSeq | 1.46.0 |
| DescTools | 0.99.50 |
| cowplot | 1.1.1 |
| base64url | 1.4 |
| ars | 0.6 |
| stats | 4.3.0 |
| wheatmap | 0.2.0 |
| uuid | 1.1-1 |
| uroot | 2.1-2 |
| TMB | 1.9.6 |
| rhdf5 | 2.44.0 |
| nloptr | 2.0.3 |
| globaltest | 5.54.0 |
| dunn.test | 1.3.5 |
| DO.db | 2.9 |
| DEsingle | 1.20.0 |
| ChAMP | 2.30.0 |
| StanHeaders | 2.26.28 |
| sjlabelled | 1.2.0 |
| rio | 1.0.1 |
| preseqR | 4.0.0 |
| performance | 0.10.5 |
| paxtoolsr | 1.34.0 |
| NLP | 0.2-1 |
| network | 1.18.1 |
| minpack.lm | 1.2-4 |
| metadat | 1.2-0 |
| impute | 1.74.1 |
| IlluminaHumanMethylation27k.db | 1.4.8 |
| ggdag | 0.2.10 |
| fdapace | 0.5.9 |
| doParallel | 1.0.17 |
| DelayedMatrixStats | 1.22.6 |
| bootstrap | 2019.6 |
| affyio | 1.70.0 |
| TCGAbiolinksGUI.data | 1.20.0 |
| sys | 3.4.2 |
| siggenes | 1.74.0 |
| sets | 1.0-24 |
| rredlist | 0.7.1 |
| ps | 1.7.5 |
| pedigreemm | 0.3-3 |
| nonnest2 | 0.5-6 |
| Matrix | 1.6-1.1 |
| logging | 0.10-108 |
| lavaan | 0.6-16 |
| ggplot.multistats | 1.0.0 |
| ggplot2 | 3.4.4 |
| GenomicScores | 2.12.1 |
| fda | 6.1.4 |
| xtable | 1.8-4 |
| VariantAnnotation | 1.46.0 |
| TFBSTools | 1.38.0 |
| spatstat.linnet | 3.1-1 |
| rstanarm | 2.26.1 |
| Rsamtools | 2.16.0 |
| phyclust | 0.1-34 |
| MBESS | 4.9.2 |
| IlluminaHumanMethylation450kanno.ilmn12.hg19 | 0.6.1 |
| GenomicRanges | 1.52.1 |
| BB | 2019.10-1 |
| xtermStyle | 3.0.5 |
| widgetTools | 1.78.0 |
| tab | 5.1.1 |
| reactlog | 1.1.1 |
| pvclust | 2.2-0 |
| processx | 3.8.2 |
| nleqslv | 3.3.4 |
| Matching | 4.10-14 |
| maftools | 2.16.0 |
| linprog | 0.9-4 |
| jomo | 2.7-6 |
| iotools | 0.3-2 |
| htmlTable | 2.4.1 |
| HSMMSingleCell | 1.20.0 |
| gtable | 0.3.4 |
| gee | 4.13-25 |
| gap.datasets | 0.0.6 |
| gamlss.dist | 6.1-1 |
| fpc | 2.2-10 |
| falconx | 0.2 |
| BSgenome.Dmelanogaster.UCSC.dm2 | 1.4.0 |
| BANOVA | 1.2.1 |
| vctrs | 0.6.4 |
| truncdist | 1.0-2 |
| seqinr | 4.2-30 |
| PolynomF | 2.0-5 |
| org.Rn.eg.db | 3.17.0 |
| lmm | 1.4 |
| ggnewscale | 0.4.9 |
| DEoptim | 2.2-8 |
| CircStats | 0.2-6 |
| taxize | 0.9.100 |
| survRM2 | 1.0-4 |
| OrganismDbi | 1.42.0 |
| OpenImageR | 1.3.0 |
| opencv | 0.2.3 |
| nucleR | 2.32.0 |
| mouse4302frmavecs | 1.5.0 |
| minfi | 1.46.0 |
| loomR | 0.2.0 |
| Iso | 0.0-21 |
| hgu133a2cdf | 2.18.0 |
| goseq | 1.52.0 |
| futile.options | 1.0.1 |
| FlowSOM | 2.8.0 |
| wikitaxa | 0.4.0 |
| shinyBS | 0.61.1 |
| segmented | 1.6-4 |
| scater | 1.28.0 |
| rstpm2 | 1.6.2 |
| reshape2 | 1.4.4 |
| rentrez | 1.2.3 |
| PKI | 0.1-12 |
| pcaPP | 2.0-3 |
| paran | 1.5.2 |
| NuPoP | 2.8.1 |
| magrittr | 2.0.3 |
| kSamples | 1.2-10 |
| KEGGdzPathwaysGEO | 1.38.0 |
| gmp | 0.7-2 |
| ggsignif | 0.6.4 |
| ggjoy | 0.4.1 |
| gbm | 2.1.8.1 |
| FME | 1.3.6.3 |
| DRR | 0.0.4 |
| Category | 2.66.0 |
| ash | 1.0-15 |
| tools | 4.3.0 |
| RTCGA | 1.30.0 |
| rmeta | 3.0 |
| RLRsim | 3.1-8 |
| repr | 1.1.6 |
| gamlss | 5.4-20 |
| curry | 0.1.1 |
| cubature | 2.1.0 |
| biocViews | 1.68.2 |
| arrayQualityMetrics | 3.56.0 |
| scDD | 1.24.0 |
| RWeka | 0.4-46 |
| R.oo | 1.25.0 |
| RhpcBLASctl | 0.23-42 |
| pbmcapply | 1.5.1 |
| IlluminaHumanMethylationEPICmanifest | 0.3.0 |
| ggrastr | 1.0.2 |
| CompQuadForm | 1.4.3 |
| biomformat | 1.28.0 |
| vegan | 2.6-4 |
| tximportData | 1.28.0 |
| scrime | 1.3.5 |
| scde | 2.28.2 |
| Rvmmin | 2018-4.17.1 |
| rjags | 4-14 |
| quantsmooth | 1.66.0 |
| mvnfast | 0.2.8 |
| muhaz | 1.2.6.4 |
| janitor | 2.2.0 |
| itertools | 0.1-3 |
| IlluminaHumanMethylationEPICanno.ilm10b2.hg19 | 0.6.0 |
| ICS | 1.4-1 |
| googleAuthR | 2.0.1 |
| corrplot | 0.92 |
| circlize | 0.4.15 |
| BiocVersion | 3.17.1 |
| visNetwork | 2.1.2 |
| tilingArray | 1.78.0 |
| swagger | 3.33.1 |
| spatstat | 3.0-6 |
| scatterpie | 0.2.1 |
| RItools | 0.3-3 |
| party | 1.3-13 |
| MendelianRandomization | 0.9.0 |
| LogicReg | 1.6.6 |
| gageData | 2.38.0 |
| earth | 5.3.2 |
| cytomapper | 1.12.0 |
| Canopy | 1.3.0 |
| urca | 1.3-3 |
| txtplot | 1.0-4 |
| SparseM | 1.81 |
| rnoaa | 1.4.0 |
| qlcMatrix | 0.9.7 |
| NBPSeq | 0.3.1 |
| msir | 1.3.3 |
| mosaicCore | 0.9.2.1 |
| iClusterPlus | 1.36.1 |
| GenomeInfoDb | 1.36.4 |
| DOT | 0.1 |
| cicero | 1.18.0 |
| trust | 0.1-8 |
| tictoc | 1.2 |
| spatstat.sparse | 3.0-2 |
| servr | 0.27 |
| Rhdf5lib | 1.22.1 |
| polspline | 1.1.23 |
| pkgdown | 2.0.7 |
| pals | 1.8 |
| later | 1.3.1 |
| kpmt | 0.1.0 |
| GEOmetadb | 1.62.0 |
| gaussquad | 1.0-3 |
| ensembldb | 2.24.1 |
| copula | 1.1-2 |
| canine2.db | 3.13.0 |
| BH | 1.81.0-1 |
| rgl | 1.2.1 |
| RANN | 2.6.1 |
| PROcess | 1.76.0 |
| geneplotter | 1.78.0 |
| fishplot | 0.5.1 |
| egg | 0.4.5 |
| ActCR | 0.3.0 |
| urlchecker | 1.0.1 |
| synchronicity | 1.3.5 |
| RcppParallel | 5.1.7 |
| qrng | 0.0-9 |
| proxy | 0.4-27 |
| pd.mouse430.2 | 3.12.0 |
| parallelly | 1.36.0 |
| metafor | 4.4-0 |
| leidenbase | 0.1.25 |
| KEGGREST | 1.40.1 |
| hapsim | 0.31 |
| GenomicDataCommons | 1.24.3 |
| cvAUC | 1.1.4 |
| basilisk | 1.12.1 |
| aod | 1.3.2 |
| zeallot | 0.1.0 |
| rsconnect | 1.1.1 |
| R.devices | 2.17.1 |
| progressr | 0.14.0 |
| patchwork | 1.1.3 |
| org.Cf.eg.db | 3.17.0 |
| inflection | 1.3.6 |
| IlluminaHumanMethylation27kmanifest | 0.4.0 |
| hgu133plus2.db | 3.13.0 |
| fds | 1.8 |
| dplyr | 1.1.3 |
| cometExactTest | 0.1.5 |
| combinat | 0.0-8 |
| cmm | 1.0 |
| methods | 4.3.0 |
| useful | 1.2.6 |
| trimcluster | 0.1-5 |
| survey | 4.2-1 |
| SQUAREM | 2021.1 |
| reticulate | 1.34.0 |
| phia | 0.2-1 |
| pfamAnalyzeR | 1.0.1 |
| ParamHelpers | 1.14.1 |
| gaston | 1.5.9 |
| clock | 0.7.0 |
| affxparser | 1.72.0 |
| weights | 1.0.4 |
| tseries | 0.10-54 |
| TFisher | 0.2.0 |
| tensorflow | 2.14.0 |
| tensor | 1.5 |
| Repitools | 1.46.0 |
| R2HTML | 2.3.3 |
| JM | 1.5-2 |
| hgu95av2probe | 2.18.0 |
| fastICA | 1.2-3 |
| compute.es | 0.2-5 |
| bigmemory | 4.6.1 |
| sourcetools | 0.1.7-1 |
| PFAM.db | 3.17.0 |
| PAIRADISE | 1.16.0 |
| lsa | 0.73.3 |
| iCNV | 1.20.0 |
| future.batchtools | 0.12.0 |
| EBImage | 4.42.0 |
| ddalpha | 1.3.13 |
| ctc | 1.74.0 |
| clusterGeneration | 1.3.8 |
| broom.mixed | 0.2.9.4 |
| AnnotationHub | 3.8.0 |
| svUnit | 1.0.6 |
| stringfish | 0.15.8 |
| shinytest | 1.5.3 |
| scattermore | 1.2 |
| mixOmics | 6.24.0 |
| misc3d | 0.9-1 |
| metap | 1.9 |
| IlluminaHumanMethylationEPICanno.ilm10b4.hg19 | 0.6.0 |
| hwriter | 1.3.2.1 |
| gmm | 1.8 |
| Glimma | 2.10.0 |
| ddCt | 1.56.0 |
| bc3net | 1.0.4 |
| WikipediR | 1.5.0 |
| RSpectra | 0.16-1 |
| RColorBrewer | 1.1-3 |
| perry | 0.3.1 |
| networkLite | 1.0.5 |
| mvord | 1.1.1 |
| minqa | 1.2.6 |
| maxLik | 1.5-2 |
| htmlwidgets | 1.6.2 |
| hgu133a2.db | 3.13.0 |
| HDF5Array | 1.28.1 |
| ggridges | 0.5.4 |
| gcrma | 2.72.0 |
| epiDisplay | 3.5.0.2 |
| bumphunter | 1.42.0 |
| compiler | 4.3.0 |
| TxDb.Hsapiens.UCSC.hg38.knownGene | 3.17.0 |
| timeDate | 4022.108 |
| sysfonts | 0.8.8 |
| statnet | 2019.6 |
| stargazer | 5.2.3 |
| S4Arrays | 1.0.6 |
| GOSim | 1.38.0 |
| getPass | 0.2-2 |
| Ecfun | 0.3-2 |
| Ecdat | 0.4-2 |
| dotCall64 | 1.1-0 |
| beanplot | 1.3.1 |
| foreign | 0.8-84 |
| roxygen2 | 7.2.3 |
| robCompositions | 2.4.1 |
| NMOF | 2.8-0 |
| limma | 3.56.2 |
| hu6800cdf | 2.18.0 |
| hierfstat | 0.5-11 |
| hgu133a2probe | 2.18.0 |
| gower | 1.0.1 |
| fontLiberation | 0.1.0 |
| flexclust | 1.4-1 |
| ca | 0.71.1 |
| WES.1KG.WUGSC | 1.32.0 |
| usethis | 2.2.2 |
| textshaping | 0.3.7 |
| shiny | 1.7.5.1 |
| R2jags | 0.7-1 |
| quadprog | 1.5-8 |
| phylobase | 0.8.10 |
| microbenchmark | 1.4.10 |
| MAVE | 1.3.11 |
| graphlayouts | 1.0.1 |
| GEOquery | 2.68.0 |
| gdsfmt | 1.36.1 |
| fitdistrplus | 1.1-11 |
| ChAMPdata | 2.32.0 |
| argparser | 0.7.1 |
| adegenet | 2.1.10 |
| threg | 1.0.3 |
| targets | 1.3.2 |
| scatterplot3d | 0.3-44 |
| promises | 1.2.1 |
| prettydoc | 0.4.1 |
| poilog | 0.4.2 |
| pkgmaker | 0.32.10 |
| pd.mogene.2.0.st | 3.14.1 |
| pbapply | 1.7-2 |
| llogistic | 1.0.3 |
| humanomni5quadv1bCrlmm | 1.0.0 |
| ellipse | 0.5.0 |
| doSNOW | 1.0.20 |
| data.tree | 1.0.0 |
| topGO | 2.52.0 |
| squash | 1.0.9 |
| spatstat.explore | 3.2-3 |
| s2 | 1.1.4 |
| Runuran | 0.38 |
| rngWELL | 0.10-9 |
| rmdformats | 1.0.4 |
| posterior | 1.4.1 |
| NbClust | 3.0.1 |
| ICC | 2.4.0 |
| fBasics | 4031.95 |
| extrafontdb | 1.0 |
| dfidx | 0.0-5 |
| xgboost | 1.7.5.1 |
| viridis | 0.6.4 |
| universalmotif | 1.18.1 |
| treeio | 1.24.3 |
| survivalAnalysis | 0.3.0 |
| stepPlr | 0.93 |
| pillar | 1.9.0 |
| parmigene | 1.1.0 |
| packrat | 0.9.2 |
| MKmisc | 1.9 |
| insight | 0.19.6 |
| hms | 1.1.3 |
| hgu95av2cdf | 2.18.0 |
| fftwtools | 0.9-11 |
| etm | 1.1.1 |
| cvTools | 0.3.2 |
| coxme | 2.2-18.1 |
| BWStest | 0.2.3 |
| bkmr | 0.2.2 |
| bindrcpp | 0.2.2 |
| aroma.light | 3.30.0 |
| XVector | 0.40.0 |
| units | 0.8-2 |
| samr | 3.0 |
| registry | 0.5-1 |
| recipes | 1.0.8 |
| intansv | 1.40.0 |
| ICSNP | 1.1-2 |
| feather | 0.3.5 |
| extrafont | 0.19 |
| base64 | 2.0.1 |
| worrms | 0.4.3 |
| wk | 0.8.0 |
| venneuler | 1.1-3 |
| VanillaICE | 1.62.0 |
| ucminf | 1.2.0 |
| svglite | 2.1.2 |
| SingleR | 2.2.0 |
| Rsolnp | 1.16 |
| reshape | 0.8.9 |
| phytools | 1.9-16 |
| mzID | 1.38.0 |
| km.ci | 0.5-6 |
| isobar | 1.46.0 |
| highr | 0.10 |
| filehash | 2.4-5 |
| densvis | 1.10.3 |
| CVST | 0.2-3 |
| arm | 1.13-1 |
| nnet | 7.3-19 |
| venn | 1.11 |
| sendmailR | 1.4-0 |
| rappdirs | 0.3.3 |
| plot3Drgl | 1.0.4 |
| pepr | 0.4.0 |
| multcomp | 1.4-25 |
| mgcv | 1.9-0 |
| ggalt | 0.4.0 |
| dtw | 1.23-1 |
| cosinor | 1.2.3 |
| broom | 1.0.5 |
| translations | 4.3.0 |
| rbibutils | 2.2.15 |
| polynom | 1.4-1 |
| org.Mm.eg.db | 3.17.0 |
| msigdbr | 7.5.1 |
| chromVAR | 1.22.1 |
| BSgenome.Hsapiens.UCSC.hg19 | 1.4.3 |
| BPSC | 0.99.2 |
| sensitivity | 1.29.0 |
| pbatR | 2.2-16 |
| mclogit | 0.9.6 |
| lmom | 3.0 |
| GlobalOptions | 0.1.2 |
| ggfun | 0.1.3 |
| geiger | 2.0.11 |
| episensr | 1.3.0 |
| compositions | 2.0-6 |
| colorspace | 2.1-0 |
| BSgenome.Scerevisiae.UCSC.sacCer1 | 1.4.0 |
| BSgenome.Hsapiens.UCSC.hg18 | 1.3.1000 |
| bold | 1.3.0 |
| splines | 4.3.0 |
| yulab.utils | 0.1.0 |
| scuttle | 1.10.3 |
| rae230aprobe | 2.18.0 |
| pvca | 1.40.0 |
| proto | 1.0.0 |
| oompaData | 3.1.3 |
| jsonlite | 1.8.7 |
| gert | 2.0.0 |
| foreach | 1.5.2 |
| forcats | 1.0.0 |
| fgsea | 1.26.0 |
| farver | 2.1.1 |
| FactoMineR | 2.9 |
| expm | 0.999-7 |
| credentials | 2.0.1 |
| BSgenome.Scerevisiae.UCSC.sacCer2 | 1.4.0 |
| BiocManager | 1.30.22 |
| subSeq | 1.30.0 |
| snowFT | 1.6-1 |
| pipeFrame | 1.16.0 |
| NMF | 0.26 |
| logicFS | 2.20.0 |
| LGEWIS | 1.1 |
| IsoformSwitchAnalyzeR | 2.0.1 |
| fit.models | 0.64 |
| ExomeDepth | 1.1.16 |
| BiocGenerics | 0.46.0 |
| codetools | 0.2-19 |
| sna | 2.7-1 |
| sgeostat | 1.0-27 |
| R.methodsS3 | 1.8.2 |
| multicool | 0.1-12 |
| MSQC | 1.1.0 |
| lmerTest | 3.1-3 |
| lazyeval | 0.2.2 |
| ergm.multi | 0.2.0 |
| ROpenCVLite | 4.70.0 |
| rBayesianOptimization | 1.2.0 |
| quantreg | 5.97 |
| plumber | 1.2.1 |
| pbs | 1.1 |
| natserv | 1.0.0 |
| naivebayes | 0.9.7 |
| monocle | 2.28.0 |
| maps | 3.4.1 |
| lmtest | 0.9-40 |
| interactiveDisplayBase | 1.38.0 |
| flexsurv | 2.2.2 |
| flexmix | 2.3-19 |
| CODEX2 | 1.3.0 |
| BSgenome.Mmusculus.UCSC.mm10 | 1.4.3 |
| bayesplot | 1.10.0 |
| accelerometry | 3.1.2 |
| WikidataR | 2.3.3 |
| TxDb.Rnorvegicus.UCSC.rn4.ensGene | 3.2.2 |
| TrajectoryUtils | 1.8.0 |
| strex | 1.6.0 |
| spls | 2.2-3 |
| shinyFiles | 0.9.3 |
| sesame | 1.18.4 |
| RInside | 0.2.18 |
| readr | 2.1.4 |
| rainbow | 3.7 |
| PMA | 1.2-2 |
| packcircles | 0.3.6 |
| IRkernel | 1.3.2 |
| IlluminaHumanMethylation27kanno.ilmn12.hg19 | 0.6.0 |
| biclust | 2.0.3.1 |
| stats4 | 4.3.0 |
| base | 4.3.0 |
| tradeSeq | 1.14.0 |
| tclust | 1.5-4 |
| synthpop | 1.8-0 |
| S4Vectors | 0.38.2 |
| rstatix | 0.7.2 |
| proj4 | 1.0-12 |
| oz | 1.0-22 |
| ordinal | 2022.11-16 |
| mitools | 2.4 |
| meta | 6.5-0 |
| googledrive | 2.1.1 |
| ggpubr | 0.6.0 |
| bitops | 1.0-7 |
| affycomp | 1.76.1 |
| WriteXLS | 6.4.0 |
| waveslim | 1.8.4 |
| TxDb.Hsapiens.UCSC.hg19.knownGene | 3.2.2 |
| tsna | 0.3.5 |
| thgenetics | 0.4-2 |
| sparseMatrixStats | 1.12.2 |
| readxl | 1.4.3 |
| RcppRoll | 0.3.0 |
| RcppProgress | 0.4.2 |
| postlogic | 0.1.0.1 |
| ggeffects | 1.3.2 |
| emmeans | 1.8.9 |
| chk | 0.9.1 |
| betareg | 3.1-4 |
| bbmle | 1.0.25 |
| BatchJobs | 1.9 |
| UpSetR | 1.4.0 |
| timsac | 1.3.8-4 |
| tergm | 4.2.0 |
| subplex | 1.8 |
| spatstat.model | 3.2-6 |
| sparsesvd | 0.2-2 |
| rversions | 2.1.2 |
| pd.genomewidesnp.6 | 3.14.1 |
| oompaBase | 3.2.9 |
| iterpc | 0.4.2 |
| irr | 0.84.1 |
| IRdisplay | 1.1 |
| gmodels | 2.18.1.1 |
| Exact | 3.2 |
| withr | 2.5.1 |
| tidyr | 1.3.0 |
| spatstat.data | 3.0-1 |
| showimage | 1.0.0 |
| RWiener | 1.3-3 |
| RSQLite | 2.3.1 |
| RcppAnnoy | 0.0.21 |
| perm | 1.0-0.4 |
| lava | 1.7.2.1 |
| intervals | 0.15.4 |
| FlowSorted.Blood.EPIC | 2.4.2 |
| EBSeq | 1.40.0 |
| zebrafishcdf | 2.18.0 |
| vroom | 1.6.4 |
| TxDb.Mmusculus.UCSC.mm9.knownGene | 3.2.2 |
| robustHD | 0.8.0 |
| RNOmni | 1.0.1.2 |
| quantmod | 0.4.25 |
| pastecs | 1.3.21 |
| nplplot | 4.6 |
| mutoss | 0.1-13 |
| hdf5r | 1.3.8 |
| glmpath | 0.98 |
| extRemes | 2.1-3 |
| chipseq | 1.50.0 |
| tsne | 0.1-3.1 |
| tripack | 1.3-9.1 |
| TailRank | 3.2.2 |
| rsample | 1.2.0 |
| qqconf | 1.3.2 |
| orthogonalsplinebasis | 0.1.7 |
| modeltools | 0.2-23 |
| maxstat | 0.7-25 |
| glmmML | 1.1.5 |
| GenomicAlignments | 1.36.0 |
| fail | 1.3 |
| enpls | 6.1 |
| brew | 1.0-8 |
| vsn | 3.68.0 |
| survminer | 0.4.9 |
| relimp | 1.0-5 |
| RcppNumerical | 0.6-0 |
| pan | 1.9 |
| mmap | 0.6-21 |
| listenv | 0.9.0 |
| gh | 1.4.0 |
| DynDoc | 1.78.0 |
| denstrip | 1.5.4 |
| CpGassoc | 2.60 |
| ConsensusClusterPlus | 1.64.0 |
| bios2mds | 1.2.3 |
| beachmat | 2.16.0 |
| tcltk | 4.3.0 |
| wdm | 0.2.4 |
| rhdf5filters | 1.12.1 |
| prabclus | 2.3-2 |
| orthopolynom | 1.0-6.1 |
| NADA | 1.6-1.1 |
| KEGGgraph | 1.60.0 |
| ggstance | 0.3.6 |
| gdtools | 0.3.4 |
| gamm4 | 0.2-6 |
| epiR | 2.0.65 |
| dygraphs | 1.1.1.6 |
| dmrseq | 1.20.1 |
| tidyverse | 2.0.0 |
| systemPipeR | 2.6.3 |
| ruv | 0.9.7.1 |
| robust | 0.7-2 |
| Rfit | 0.24.6 |
| RefPlus | 1.70.0 |
| mutSignatures | 2.1.1 |
| minet | 3.58.0 |
| hgu95av2 | 2.2.0 |
| getopt | 1.20.4 |
| geomorph | 4.0.6 |
| fftw | 1.0-7 |
| falcon | 0.2 |
| evd | 2.3-6.1 |
| diptest | 0.76-0 |
| dimRed | 0.2.6 |
| depmap | 1.14.0 |
| dendextend | 1.17.1 |
| conditionz | 0.1.0 |
| ChromHeatMap | 1.54.0 |
| chopsticks | 1.66.0 |
| arrangements | 1.1.9 |
| truncnorm | 1.0-9 |
| systemfit | 1.1-30 |
| RApiSerialize | 0.1.2 |
| randomForestSRC | 3.2.2 |
| plot3D | 1.4 |
| optparse | 1.7.3 |
| optimParallel | 1.0-2 |
| mclust | 6.0.0 |
| manipulateWidget | 0.11.1 |
| logistf | 1.26.0 |
| hash | 2.2.6.3 |
| ggrepel | 0.9.4 |
| ggbio | 1.48.0 |
| EMCluster | 0.2-15 |
| coxphf | 1.13.4 |
| ClusterR | 1.3.1 |
| brms | 2.20.4 |
| bridgesampling | 1.1-2 |
| ASSET | 2.18.0 |
| webutils | 1.1 |
| uwot | 0.1.16 |
| SKAT | 2.2.5 |
| shinystan | 2.6.0 |
| scran | 1.28.2 |
| mnormt | 2.1.1 |
| kernlab | 0.9-32 |
| HIBAG | 1.36.4 |
| h2o | 3.42.0.2 |
| grpreg | 3.4.0 |
| gbRd | 0.4-11 |
| clisymbols | 1.2.0 |
| supraHex | 1.38.0 |
| Seurat | 4.4.0 |
| sctransform | 0.4.1 |
| rtracklayer | 1.60.1 |
| Rmisc | 1.5.1 |
| reprex | 2.0.2 |
| RcisTarget | 1.20.0 |
| officer | 0.6.2 |
| micEcon | 0.6-18 |
| callr | 3.7.3 |
| XML | 3.99-0.14 |
| Rttf2pt1 | 1.3.12 |
| restfulr | 0.0.15 |
| ratelimitr | 0.4.1 |
| nnls | 1.5 |
| mlr | 2.19.1 |
| maptools | 1.1-8 |
| esATAC | 1.22.0 |
| crlmm | 1.58.0 |
| cobs | 1.3-5 |
| splines2 | 0.5.1 |
| RJSONIO | 1.3-1.8 |
| RcmdrPlugin.TeachingDemos | 1.2-0 |
| RcmdrMisc | 2.9-1 |
| Rcmdr | 2.9-1 |
| oligo | 1.64.1 |
| mvbutils | 2.8.232 |
| irlba | 2.3.5.1 |
| gtools | 3.9.4 |
| Epi | 2.47.1 |
| eha | 2.11.1 |
| cubelyr | 1.0.2 |
| clipr | 0.8.0 |
| cghMCR | 1.58.0 |
| Biostrings | 2.68.1 |
| aplot | 0.2.2 |
| zip | 2.3.0 |
| yaml | 2.3.7 |
| xlsx | 0.6.5 |
| Rdpack | 2.5 |
| rae230acdf | 2.18.0 |
| networkD3 | 0.4 |
| mets | 1.3.2 |
| genefilter | 1.82.1 |
| formula.tools | 1.7.1 |
| epitools | 0.5-10.1 |
| EDASeq | 2.34.0 |
| dynamicTreeCut | 1.63-1 |
| distributional | 0.3.2 |
| datawizard | 0.9.0 |
| cluster | 2.1.4 |
| signatureSearch | 1.14.0 |
| sfsmisc | 1.1-16 |
| seewave | 2.2.3 |
| Rhtslib | 2.2.0 |
| RcppZiggurat | 0.1.6 |
| qap | 0.1-2 |
| pander | 0.6.5 |
| optextras | 2019-12.4 |
| ini | 0.3.1 |
| Illumina450ProbeVariants.db | 1.36.0 |
| frma | 1.52.0 |
| ExperimentHub | 2.8.1 |
| EnvStats | 2.8.1 |
| dtwclust | 5.5.12 |
| doFuture | 1.0.0 |
| SuppDists | 1.1-9.7 |
| supclust | 1.1-1 |
| slam | 0.1-50 |
| shinythemes | 1.2.0 |
| SeuratObject | 4.1.4 |
| rprojroot | 2.0.3 |
| PSCBS | 0.66.0 |
| pixmap | 0.4-12 |
| nycflights13 | 1.0.2 |
| HardyWeinberg | 1.7.5 |
| ggtext | 0.1.2 |
| ecodist | 2.0.9 |
| EBarrays | 2.64.0 |
| devEMF | 4.4-1 |
| cOde | 1.1.1 |
| cachem | 1.0.8 |
| xmlparsedata | 1.0.5 |
| TFMPvalue | 0.0.9 |
| ROC | 1.76.0 |
| princurve | 2.1.6 |
| hoardr | 0.5.3 |
| gRbase | 2.0.0 |
| GeneRegionScan | 1.56.0 |
| fmsb | 0.7.5 |
| fastmap | 1.1.1 |
| fastcluster | 1.2.3 |
| DT | 0.30 |
| doBy | 4.6.19 |
| verification | 1.42 |
| SpatialExperiment | 1.10.0 |
| ROpenCV | unknown |
| rootSolve | 1.8.2.4 |
| RcppHNSW | 0.5.0 |
| qtl | 1.60 |
| kde1d | 1.0.5 |
| httpuv | 1.6.11 |
| HMMcopy | 1.42.0 |
| fission | 1.20.0 |
| fibroEset | 1.42.0 |
| emdbook | 1.3.13 |
| drc | 3.0-1 |
| broman | 0.80 |
| Biobase | 2.60.0 |
| bigmemory.sri | 0.1.6 |
| assertthat | 0.2.1 |
| V8 | 4.4.0 |
| udunits2 | 0.13.2.1 |
| selectr | 0.4-2 |
| ROCR | 1.0-11 |
| permute | 0.9-7 |
| pcaMethods | 1.92.0 |
| np | 0.60-17 |
| mvmeta | 1.0.3 |
| interp | 1.1-4 |
| hthgu133aprobe | 2.18.0 |
| GenSA | 1.1.10.1 |
| FDb.InfiniumMethylation.hg19 | 2.2.0 |
| fastDummies | 1.7.3 |
| dvmisc | 1.1.4 |
| directlabels | 2023.8.25 |
| car | 3.1-2 |
| BiocNeighbors | 1.18.0 |
| AnnotationDbi | 1.62.2 |
| animation | 2.7 |
| wordcloud | 2.6 |
| sandwich | 3.0-2 |
| pspline | 1.0-19 |
| methylKit | 1.26.0 |
| memoise | 2.0.1 |
| IlluminaHumanMethylation450kmanifest | 0.4.0 |
| FlowSorted.CordBlood.450k | 1.28.0 |
| FDb.InfiniumMethylation.hg18 | 2.2.0 |
| DOSE | 3.26.2 |
| debugme | 1.1.0 |
| corrgram | 1.14 |
| c3net | 1.1.1.1 |
| zoo | 1.8-12 |
| triangle | 1.0 |
| tis | 1.39 |
| splancs | 2.01-44 |
| SimpleITK | 2.1.1.1 |
| scRNAseq | 2.14.0 |
| pd.mogene.1.0.st.v1 | 3.14.1 |
| parameters | 0.21.2 |
| mcmcplots | 0.4.3 |
| graphite | 1.46.0 |
| GGally | 2.1.2 |
| actuar | 3.3-2 |
| utf8 | 1.2.3 |
| TxDb.Hsapiens.UCSC.hg18.knownGene | 3.2.2 |
| simsurv | 1.0.0 |
| rncl | 0.8.7 |
| RcppThread | 2.1.6 |
| poweRlaw | 0.70.6 |
| GSVA | 1.48.3 |
| flowCore | 2.12.2 |
| ergm | 4.5.0 |
| dtplyr | 1.3.1 |
| collapse | 2.0.3 |
| bezier | 1.1.2 |
| bayesm | 3.1-6 |
| AnnotationFilter | 1.24.0 |
| svgPanZoom | 0.3.4 |
| ResourceSelection | 0.3-6 |
| philentropy | 0.7.0 |
| mice | 3.16.0 |
| lubridate | 1.9.3 |
| lme4 | 1.1-34 |
| hardhat | 1.3.0 |
| GO.db | 3.17.0 |
| GMMAT | 1.4.1 |
| geepack | 1.3.9 |
| Formula | 1.2-5 |
| DSS | 2.48.0 |
| cummeRbund | 2.42.0 |
| catmap | 1.6.4 |
| splitstackshape | 1.4.8 |
| RnBeads | 2.18.1 |
| RcppEigen | 0.3.3.9.3 |
| modeldata | 1.2.0 |
| MLInterfaces | 1.80.0 |
| ISwR | 2.0-8 |
| ineq | 0.2-13 |
| glue | 1.6.2 |
| GenomicFeatures | 1.52.2 |
| doMC | 1.3.8 |
| dendsort | 0.3.4 |
| dbplyr | 2.3.4 |
| checkmate | 2.2.0 |
| calibrate | 1.7.7 |
| BSgenome.Mmusculus.UCSC.mm9 | 1.4.0 |
| tidytree | 0.4.5 |
| stringi | 1.7.12 |
| rARPACK | 0.11-0 |
| plogr | 0.2.0 |
| nabor | 0.5.0 |
| klaR | 1.7-2 |
| infotheo | 1.2.0.1 |
| glmmTMB | 1.1.8 |
| ggmap | 3.0.2 |
| genalg | 0.2.1 |
| effectsize | 0.8.6 |
| dynlm | 0.3-6 |
| Deriv | 4.1.3 |
| cytolib | 2.12.1 |
| crul | 1.4.0 |
| AlgDesign | 1.2.1 |
| TxDb.Dmelanogaster.UCSC.dm3.ensGene | 3.2.2 |
| tmvtnorm | 1.5 |
| tfautograph | 0.3.2 |
| sodium | 1.2.1 |
| regioneR | 1.32.0 |
| PureCN | 2.6.4 |
| pbivnorm | 0.6.0 |
| networkDynamic | 0.11.3 |
| leaps | 3.1 |
| gistr | 0.9.0 |
| gam | 1.22-2 |
| fs | 1.6.3 |
| ffpe | 1.44.0 |
| fastseg | 1.46.0 |
| DirichletMultinomial | 1.42.0 |
| agricolae | 1.3-6 |
| spatstat.geom | 3.2-7 |
| sass | 0.4.7 |
| RMariaDB | 1.3.0 |
| read.gt3x | 1.2.0 |
| pscl | 1.5.5.1 |
| numbers | 0.8-5 |
| multtest | 2.56.0 |
| isoband | 0.2.7 |
| GIGrvg | 0.8 |
| genoPlotR | 0.8.11 |
| effects | 4.2-2 |
| distr | 2.9.2 |
| devtools | 2.4.5 |
| deSolve | 1.38 |
| crayon | 1.5.2 |
| ChIPQC | 1.36.1 |
| biomaRt | 2.56.1 |
| Wrench | 1.18.0 |
| timechange | 0.2.0 |
| stabs | 0.6-4 |
| ReportingTools | 2.39.0 |
| qcc | 2.7 |
| preprocessCore | 1.62.1 |
| mixtools | 2.0.0 |
| lumi | 2.52.0 |
| igraph | 1.5.1 |
| GSEABase | 1.62.0 |
| GSA | 1.03.2 |
| GreyListChIP | 1.32.1 |
| glasso | 1.11 |
| GenomicDistributions | 1.8.0 |
| conquer | 1.3.3 |
| classInt | 0.4-10 |
| ChIPseqR | 1.54.0 |
| cellranger | 1.1.0 |
| arrow | 12.0.1 |
| apcluster | 1.4.11 |
| admisc | 0.33 |
| survMisc | 0.5.6 |
| SuperLearner | 2.0-28.1 |
| snakecase | 0.11.1 |
| sm | 2.2-5.7.1 |
| SeqArray | 1.40.1 |
| RSNNS | 0.4-16 |
| metapod | 1.8.0 |
| manipulate | 1.0.1 |
| goftest | 1.2-3 |
| diffobj | 0.3.5 |
| CNTools | 1.56.0 |
| BSgenome.Hsapiens.1000genomes.hs37d5 | 0.99.1 |
| BRAIN | 1.46.0 |
| BBmisc | 1.13 |
| AllelicImbalance | 1.38.0 |
| ScaledMatrix | 1.8.1 |
| ritis | 1.0.0 |
| RcppGSL | 0.3.13 |
| parallelMap | 1.5.1 |
| mouse4302cdf | 2.18.0 |
| merTools | 0.6.1 |
| ieugwasr | 0.1.5 |
| HSAUR2 | 1.1-20 |
| GenomeInfoDbData | 1.2.10 |
| fontquiver | 0.2.1 |
| EGSEAdata | 1.28.0 |
| dbscan | 1.1-11 |
| caret | 6.0-94 |
| BiocBaseUtils | 1.2.0 |
| BiasedUrn | 2.0.11 |
| tmvnsim | 1.0-2 |
| snowfall | 1.84-6.2 |
| signatureSearchData | 1.14.0 |
| rematch | 2.0.0 |
| Rbowtie2 | 2.6.0 |
| randomForest | 4.7-1.1 |
| MLmetrics | 1.1.1 |
| labeling | 0.4.3 |
| fastGHQuad | 1.0.1 |
| BubbleTree | 2.30.0 |
| bdsmatrix | 1.3-6 |
| afex | 1.3-0 |
| xts | 0.13.1 |
| tidytidbits | 0.3.2 |
| sn | 2.1.1 |
| rstantools | 2.3.1.1 |
| ranger | 0.15.1 |
| operator.tools | 1.6.3 |
| OpenMx | 2.21.8 |
| motifStack | 1.44.1 |
| lpsymphony | 1.28.1 |
| forecast | 8.21.1 |
| diagram | 1.6.5 |
| nlme | 3.1-162 |
| statnet.common | 4.9.0 |
| sjstats | 0.18.2 |
| rvcheck | 0.2.1 |
| RProtoBufLib | 2.12.1 |
| poorman | 0.2.6 |
| mda | 0.5-4 |
| matrixStats | 1.0.0 |
| hgu133plus2cdf | 2.18.0 |
| girafe | 1.52.0 |
| ggvis | 0.4.8 |
| BSgenome.Celegans.UCSC.ce2 | 1.4.0 |
| beeswarm | 0.4.0 |
| annotate | 1.78.0 |
| affycoretools | 1.72.0 |
| singscore | 1.20.0 |
| sem | 3.1-15 |
| pseval | 1.3.1 |
| mitml | 0.4-5 |
| JGR | 1.9-2 |
| heatmaply | 1.5.0 |
| argparse | 2.2.2 |
| VennDiagram | 1.7.3 |
| TTR | 0.24.3 |
| sjPlot | 2.8.15 |
| runjags | 2.2.1-7 |
| Rmpi | 0.7-1 |
| kinship2 | 1.9.6 |
| here | 1.0.1 |
| gpls | 1.72.0 |
| GENESIS | 2.30.0 |
| forestplot | 3.1.3 |
| doRNG | 1.8.6 |
| chron | 2.3-61 |
| arrayQuality | 1.78.0 |
| xlsxjars | 0.6.1 |
| qgraph | 1.9.5 |
| PearsonDS | 1.3.0 |
| MCMCpack | 1.6-3 |
| hmmm | 1.0-4 |
| haplo.stats | 1.9.3 |
| gitcreds | 0.1.2 |
| desc | 1.4.2 |
| BSgenome.Hsapiens.UCSC.hg38 | 1.4.5 |
| Brobdingnag | 1.2-9 |
| BiocFileCache | 2.8.0 |
| babelgene | 22.9 |
| sitmo | 2.0.2 |
| RUVSeq | 1.34.0 |
| rsvg | 2.6.0 |
| pracma | 2.4.2 |
| pacman | 0.5.1 |
| mogene20sttranscriptcluster.db | 8.8.0 |
| MESS | 0.5.12 |
| grImport | 0.9-7 |
| commonsMath | 1.2.8 |
| blme | 1.0-5 |
| askpass | 1.2.0 |
| TSP | 1.2-4 |
| SnowballC | 0.7.1 |
| SC3 | 1.28.3 |
| satuRn | 1.8.0 |
| RBGL | 1.76.0 |
| R2OpenBUGS | 3.2-3.2.1 |
| purrr | 1.0.2 |
| lpSolve | 5.6.19 |
| lambda.r | 1.2.4 |
| labelled | 2.12.0 |
| influenceR | 0.1.5 |
| ggsci | 3.0.0 |
| doMPI | 0.2.2 |
| docopt | 0.7.1 |
| dlm | 1.1-6 |
| DelayedArray | 0.26.7 |
| cyclocomp | 1.1.1 |
| clusterProfiler | 4.8.3 |
| SRAdb | 1.62.0 |
| signeR | 2.2.1 |
| schoolmath | 0.4.2 |
| SCATEData | 1.10.0 |
| rlang | 1.1.1 |
| RcppArmadillo | 0.12.6.4.0 |
| R2admb | 0.7.16.3 |
| qqman | 0.1.9 |
| mlogit | 1.1-1 |
| locfdr | 1.1-8 |
| LearnBayes | 2.15.1 |
| lawstat | 3.6 |
| kknn | 1.3.1 |
| karyoploteR | 1.26.0 |
| gclus | 1.3.2 |
| curl | 5.1.0 |
| covr | 3.6.3 |
| cancerTiming | 3.1.8 |
| aplpack | 1.3.5 |
| acepack | 1.4.2 |
| solrium | 1.2.0 |
| sf | 1.0-14 |
| rsvd | 1.0.5 |
| rmutil | 1.1.10 |
| rlemon | 0.2.1 |
| mcbiopi | 1.1.6 |
| logitnorm | 0.8.38 |
| lintr | 3.1.0 |
| leiden | 0.4.3 |
| kableExtra | 1.3.4 |
| GWASTools | 1.46.0 |
| GGIR | 3.0-0 |
| BradleyTerry2 | 1.1-2 |
| BeadDataPackR | 1.52.0 |
| bamsignals | 1.32.0 |
| Amelia | 1.8.1 |
| grid | 4.3.0 |
| shinyjs | 2.1.0 |
| sessioninfo | 1.2.2 |
| multcompView | 0.1-9 |
| mime | 0.12 |
| MALDIquant | 1.22.1 |
| Hmisc | 5.1-1 |
| CBPS | 0.23 |
| BMA | 3.18.17 |
| BiocStyle | 2.28.1 |
| AER | 1.2-10 |
| WikidataQueryServiceR | 1.0.0 |
| tm | 0.7-11 |
| SomaticSignatures | 2.36.0 |
| scales | 1.2.1 |
| ragg | 1.2.6 |
| nortest | 1.0-4 |
| MRInstruments | 0.3.2 |
| metagenomeSeq | 1.42.0 |
| isdparser | 0.4.0 |
| gplots | 3.1.3 |
| formatR | 1.14 |
| evaluate | 0.22 |
| ellipsis | 0.3.2 |
| amap | 0.8-19 |
| WGCNA | 1.72-1 |
| superpc | 1.12 |
| startupmsg | 0.9.6 |
| spelling | 2.2.1 |
| shinycssloaders | 1.0.0 |
| seqbias | 1.48.0 |
| rms | 6.7-1 |
| RcppML | 0.3.7 |
| plotrix | 3.8-2 |
| MVA | 1.0-8 |
| memisc | 0.99.31.6 |
| LaplacesDemon | 16.1.6 |
| Homo.sapiens | 1.3.1 |
| HKprocess | 0.1-1 |
| hapmapsnp6 | 1.42.0 |
| gridExtra | 2.3 |
| GPfit | 1.0-8 |
| geeM | 0.10.1 |
| FSA | 0.9.5 |
| crosstalk | 1.2.0 |
| affy | 1.78.2 |
| waldo | 0.5.1 |
| splus2R | 1.3-3 |
| smoother | 1.1 |
| SIS | 0.8-8 |
| Signac | 1.11.0 |
| renv | 1.0.3 |
| prodlim | 2023.08.28 |
| ppclust | 1.1.0 |
| MutationalPatterns | 3.10.0 |
| matlab | 1.0.4 |
| IRanges | 2.34.1 |
| igraphdata | 1.0.1 |
| Gmisc | 3.0.3 |
| ggplotify | 0.1.2 |
| geneLenDataBase | 1.36.0 |
| etrunct | 0.1 |
| entropy | 1.3.1 |
| coda | 0.19-4 |
| Affymoe4302Expr | 1.38.0 |
| tfruns | 1.5.1 |
| rstudioapi | 0.15.0 |
| RiboProfiling | 1.30.0 |
| rgeos | 0.6-2 |
| pkgconfig | 2.0.3 |
| pedgene | 3.6 |
| nor1mix | 1.3-0 |
| JASPAR2016 | 1.28.0 |
| htmltools | 0.5.6.1 |
| globals | 0.16.2 |
| ggExtra | 0.10.1 |
| ff | 4.0.9 |
| abind | 1.4-5 |
| rpart | 4.1.19 |
| umap | 0.2.10.0 |
| R.filesets | 2.15.0 |
| Nozzle.R1 | 1.1-1.1 |
| mvtnorm | 1.2-3 |
| MRMix | 0.1.0 |
| MatrixGenerics | 1.12.3 |
| googleVis | 0.7.1 |
| future | 1.33.0 |
| EnsDb.Hsapiens.v86 | 2.99.0 |
| DiagrammeR | 1.0.10 |
| convert | 1.76.0 |
| lattice | 0.21-8 |
| XLConnect | 1.0.7 |
| tkrplot | 0.0-27 |
| sesameData | 1.18.0 |
| seqminer | 9.1 |
| Rtsne | 0.16 |
| pso | 1.0.4 |
| logisticPCA | 0.2 |
| GWASdata | 1.38.1 |
| ggseqlogo | 0.1 |
| FD | 1.0-12.1 |
| BSgenome | 1.68.0 |
| arrayhelpers | 1.1-0 |
| tfestimators | 1.9.2 |
| sjmisc | 2.8.9 |
| shape | 1.4.6 |
| R.utils | 2.12.2 |
| lhs | 1.1.6 |
| httpcode | 0.3.0 |
| gridBase | 0.4-7 |
| genomeIntervals | 1.56.0 |
| furrr | 0.3.1 |
| bookdown | 0.36 |
| backports | 1.4.1 |
| affyPLM | 1.76.1 |
| tcltk2 | 1.2-11 |
| showtextdb | 3.0 |
| shinyWidgets | 0.8.0 |
| Rsubread | 2.14.2 |
| rgexf | 0.16.2 |
| rgenoud | 5.9-0.3 |
| projpred | 2.7.0 |
| polyclip | 1.10-6 |
| microbiome | 1.22.0 |
| magick | 2.7.4 |
| ltsa | 1.4.6 |
| InteractionSet | 1.28.1 |
| Icens | 1.72.0 |
| hgu133a.db | 3.13.0 |
| GENEAread | 2.0.9 |
| estimability | 1.4.1 |
| digest | 0.6.33 |
| BiocParallel | 1.34.2 |
| writexl | 1.4.2 |
| TxDb.Mmusculus.UCSC.mm10.knownGene | 3.10.0 |
| sp | 2.1-1 |
| SingleCellExperiment | 1.22.0 |
| Rgraphviz | 2.44.0 |
| ncvreg | 3.14.1 |
| mscstts | 0.6.3 |
| mixmeta | 1.2.0 |
| HDO.db | 0.99.1 |
| FateID | 0.2.2 |
| clValid | 0.7 |
| anytime | 0.3.9 |
| xaringan | 0.28 |
| webshot | 0.5.5 |
| terra | 1.7-55 |
| survivalROC | 1.0.3.1 |
| RgoogleMaps | 1.4.5.3 |
| rBiopaxParser | 2.40.0 |
| pseudo | 1.4.3 |
| plsVarSel | 0.9.10 |
| pkgcond | 0.1.1 |
| ModelMetrics | 1.2.2.2 |
| knitr | 1.44 |
| jpeg | 0.1-10 |
| JASPAR2018 | 1.1.1 |
| gridSVG | 1.7-5 |
| FNN | 1.1.3.2 |
| exactRankTests | 0.8-35 |
| changepoint | 2.2.4 |
| carData | 3.0-5 |
| biovizBase | 1.48.0 |
| xml2 | 1.3.5 |
| tximport | 1.28.0 |
| timereg | 2.0.5 |
| texreg | 1.38.6 |
| stringr | 1.5.0 |
| praznik | 11.0.0 |
| phangorn | 2.11.1 |
| OmicCircos | 1.38.0 |
| jquerylib | 0.1.4 |
| flextable | 0.9.3 |
| BSgenome.Ecoli.NCBI.20080805 | 1.3.1000 |
| bluster | 1.10.0 |
| VIM | 6.2.2 |
| tuneR | 1.4.5 |
| tibble | 3.2.1 |
| RNeXML | 2.4.11 |
| Rglpk | 0.6-5 |
| penalized | 0.9-52 |
| munsell | 0.5.0 |
| HiClimR | 2.2.1 |
| gProfileR | 0.7.0 |
| gnm | 1.1-5 |
| gap | 1.5-3 |
| facets | 0.6.2 |
| energy | 1.7-11 |
| dgof | 1.4 |
| Cairo | 1.6-1 |
| BSgenome.Cfamiliaris.UCSC.canFam2 | 1.4.0 |
| ASCAT | 3.1.2 |
| ActiveDriver | 1.0.0 |
| sloop | 1.0.1 |
| R.huge | 0.10.0 |
| Rcpp | 1.0.11 |
| png | 0.1-8 |
| missMethyl | 1.34.0 |
| MatrixModels | 0.5-2 |
| gsmoothr | 0.1.7 |
| forge | 0.2.0 |
| fastmatch | 1.1-4 |
| CATALYST | 1.24.0 |
| BiocIO | 1.10.0 |
| aroma.apd | 0.7.0 |
| aroma.affymetrix | 3.2.1 |
| aCGH | 1.78.0 |
| triebeard | 0.4.1 |
| subselect | 0.15.5 |
| SNPRelate | 1.34.1 |
| showtext | 0.9-6 |
| ShortRead | 1.58.0 |
| QuickJSR | 1.0.7 |
| profileModel | 0.6.1 |
| mixdist | 0.5-5 |
| miscF | 0.1-5 |
| mi | 1.1 |
| mboost | 2.9-8 |
| genomewidesnp6Crlmm | 1.0.7 |
| generics | 0.1.3 |
| float | 0.3-1 |
| ergm.count | 4.1.1 |
| EpiDynamics | 0.3.1 |
| enrichplot | 1.20.3 |
| Deducer | 0.7-9 |
| cleaver | 1.38.0 |
| bnclassify | 0.4.7 |
| systemfonts | 1.0.5 |
| stringdist | 0.9.10 |
| som | 0.3-5.1 |
| setRNG | 2022.4-1 |
| mlbench | 2.1-3.1 |
| lpSolveAPI | 5.5.2.0-17.10 |
| gridGraphics | 0.5-1 |
| ggforce | 0.4.1 |
| FlowSorted.Blood.450k | 1.38.0 |
| findpython | 1.0.8 |
| fields | 15.2 |
| downlit | 0.4.3 |
| aws | 2.5-3 |
| annotatr | 1.26.0 |





Manage your own packages
[top](#top)
#### Per-user R library


Users can install their own packages. By default, on RHEL8, this private library
is located at `/data/$USER/R/rhel8/%v` where `%v`
is the major.minor version of R (e.g. 4.3). This is a change from the
behaviour on RHEL7 where the default location was `~/R/%v/library`.

 Note that this directory is not created automatically by some
versions of R so it is safest to create it manually before installing R
packages.


Users can choose alternative locations for this directory by setting and
exporting the environment variable `$R_LIBS_USER` in your shell
startup script. If you are a bash user, for example, you could add the
following line to your `~/.bash_profile` to relocated your R
library:


```

export R_LIBS_USER="/data/$USER/code/R/rhel8/%v"

```

Here is an example using the [pacman](http://trinker.github.io/pacman/) 
package for easier package management:



```

[user@cn3144 ~]$ R

R version 4.1.0 (2021-05-18) -- "Camp Pontanezen"
Copyright (C) 2021 The R Foundation for Statistical Computing
Platform: x86_64-pc-linux-gnu (64-bit)

> **library(pacman)**
> **p\_isinstalled(rapport)**
[1] FALSE
> **p\_install(rapport)**
Installing package into ‘/spin1/home/linux/user/R/4.1/library’
(as ‘lib’ is unspecified)
also installing the dependency ‘rapportools’

trying URL 'http://cran.rstudio.com/src/contrib/rapportools_1.0.tar.gz'

[...snip...]
rapport installed
>

```

#### More granular (and reproducible) package management


A better approach than relying on packages installed centrally or in
your home directory is to create isolated, per project package sets. This
increases reproducibility at the cost of increased storage and potential
package installation headaches. Some packages to implement this:


* Recommended: The [renv](https://rstudio.github.io/renv/articles/renv.html) package
 is a more modern replacement for packrat.
* The older [packrat](https://rstudio.github.io/packrat/) package implements isolated per
project R libraries. See the [packrat walkthrough](https://rstudio.github.io/packrat/walkthrough.html)
for a basic introduction.



R batch job
[top](#top)
R batch jobs are similar to any other batch job. A batch script ('rjob.sh')
is created that sets up the environment and runs the R code:



```

#!/bin/bash

module load R/4.2
R --no-echo --no-restore --no-save < /data/user/Rtests/Rtest.r > /data/user/Rtests/Rtest.out

```

or use `Rscript` instead



```

#!/bin/bash

module load R/3.5
Rscript /data/user/Rtests/Rtest.r > /data/user/Rtests/Rtest.out

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] rjob.sh
```

#### Command line arguments for R scripts


R scripts can be written to accept command line arguments. The simplest
way of doing this is with the `commandArgs()` function. For
example the script 'simple\_args.R'



```

args <- commandArgs(trailingOnly=TRUE)

i <- 0
for (arg in args) {
    i <- i + 1
    cat(sprintf("arg %02i: '%s'\n", i, arg))
}

```

can be called like this



```

[user@cn3144]$ **module load R**
[user@cn3144]$ **Rscript simple.R this is a test**
arg 01: 'this'
arg 02: 'is'
arg 03: 'a'
arg 04: 'test'
[user@cn3144]$ **Rscript simple.R 'this is a test'**
arg 01: 'this is a test'
[user@cn3144]$ **R --no-echo --no-restore --no-save --args 'this is a test' < simple.R**
arg 01: 'this is a test'

```

Alternatively, commandline arguments can be parsed using the getopt package. For example:



```

library(getopt)

###
### Describe the expected command line arguments
###
# mask: 0=no argument
#       1=required argument
#       2=optional argument
spec <- matrix(c(
# long name  short name  mask  type          description(optional)
# ---------  ----------  ----  ------------  ---------------------
  'file'   , 'f',          1,  'character',  'input file',
  'verbose', 'v',          0,  'logical',    'verbose output', 
  'help'   , 'h',          0,  'logical',    'show this help message'
), byrow=TRUE, ncol=5);

# parse the command line
opt <- getopt(spec);

# show help if requested
if (!is.null(opt$help)) {
  cat(getopt(spec, usage=TRUE));
  q();
}

# set defaults
if ( is.null(opt$file) )    { opt$file    = 'testfile' }
if ( is.null(opt$verbose) ) { opt$verbose = FALSE }
print(opt)

```

This script an be used as follows



```

[user@cn3144]$ **Rscript getopt\_example.R --file some.txt --verbose**
$ARGS
character(0)

$file
[1] "some.txt"

$verbose
[1] TRUE

[user@cn3144]$ **Rscript getopt\_example.R --file some.txt**
$ARGS
character(0)

$file
[1] "some.txt"

$verbose
[1] FALSE

[user@cn3144]$ **Rscript getopt\_example.R --help**
Usage: getopt_example.R [-[-file|f] ] [-[-verbose|v]] [-[-help|h]]
 -f|--file input file
 -v|--verbose verbose output
 -h|--help show this help message

```

getopt does not have support for mixing flags and positional arguments.
There are other packages with different features and approaches that can
be used to design command line interfaces for R scripts.



Swarm of R jobs
[top](#top)
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a
set of independent commands requiring identical resources.


Create a swarmfile (e.g. rjobs.swarm). For example:



```

Rscript /data/user/R/R1  > /data/user/R/R1.out
Rscript /data/user/R/R2  > /data/user/R/R2.out
Rscript /data/user/R/R3  > /data/user/R/R3.out

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f TEMPLATE.swarm [-g #] [-t #] --module R/3.5
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g #  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t #  Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module TEMPLATE  Loads the TEMPLATE module for each subjob in the swarm 
 | |
 | |
 | |



Rswarm
[top](#top)
`Rswarm` is a utility to create a series of R input files from a
single R (master) template file with different output filenames and with unique
random number generator seeds. It will simultaneously create a swarm command
file that can be used to submit the swarm of R jobs. Rswarm was originally
developed by Lori Dodd and Trevor Reeve with modifications by the Biowulf
staff.


Say, for example, that the goal of a simulation study is to evaluate
properties of the t-test. The function "sim.fun" in file "sim.R" below
repeatedly generates random normal data with a given mean, performs a one
sample t-test (i.e. testing if the mean is different from 0), and records the
p-values.



```

#######################################
# n.samp:  size of samples generated for each simulation
# mu:      mean
# sd:      standard deviation
# nsim:    the number of simulations
# output1: output table
# seed:    the seed for set.seed
#######################################
sim.fun <- function(n.samp=100, mu=0, sd=1, n.sim, output1, seed){

    set.seed(seed)

    p.values <- c()
    for (i in 1:n.sim){
        x <- rnorm(n.samp, mean=mu, sd=sd)
        p.values <- c(p.values, t.test(x)$p.value)
    }
    saveRDS(p.values, file=output1)
}

```

To use Rswarm, create a wrapper script similar to the following ("rfile.R")



```

source("sim.R")
sim.fun(n.sim=DUMX, output1="DUMY1",seed=DUMZ)

```

using the the dummy variables which will be replaced by Rswarm.





| **Dummy variable** | **Replaced with** |
| --- | --- |
| DUMX | Number of simulations to be specified in each replicate file |
| DUMY1 | Output file 1 |
| DUMY2 | Output file 2 (optional) |
| DUMZ | Random seed |



To swarm this code, we need replicates of the rfile.R file, each with a
different seed and different output file. The Rswarm utility will create the
specified number of replicates, supply each with a different seed (from an
external file containing seed numbers), and create unique output files for each
replicate. Note, that we allow for you to specify the number of simulations
within each file, in addition to specifying the number of replicates.


For example, the following Rswarm command at the Biowulf prompt will create
2 replicate files, each specifying 50 simulations, a different seed from a
file entitled, "seedfile.txt," and unique output files.



```

[user@biowulf]$ **ls -lh**
total 8.0K
-rw-r--r-- 1 user group  63 Apr 25 12:34 rfile.R
-rw-r--r-- 1 user group 564 Apr 25 12:15 seedfile.txt
-rw-r--r-- 1 user group 547 Apr 25 12:04 sim.R
[user@biowulf]$ **head -n2 seedfile.txt**
24963
27507
[user@biowulf]$ **Rswarm --rfile=rfile.R --sfile=seedfile.txt --path=. \
 --reps=2 --sims=50 --start=0 --ext1=.rds**
The template file is rfile.R
The seed file is seedfile.txt
The path is .
The number of replicates desired is 2
The number of sims per file is 50
The starting file number is 0+1
The extension for output files 1 is .rds
The extension for output files 2 is .std.txt
Is this correct (y or n)? : y
Creating file number 1: ./rfile1.R with output ./rfile1.rds ./rfile1.std.txt and seed 24963
Creating file number 2: ./rfile2.R with output ./rfile2.rds ./rfile2.std.txt and seed 27507
[user@biowulf]$ **ls -lh**
total 16K
-rw-r--r-- 1 user group  69 Apr 25 12:39 rfile1.R
-rw-r--r-- 1 user group  69 Apr 25 12:39 rfile2.R
-rw-r--r-- 1 user group  63 Apr 25 12:34 rfile.R
-rw-r--r-- 1 user group  50 Apr 25 12:39 rfile.sw
-rw-r--r-- 1 user group 564 Apr 25 12:15 seedfile.txt
-rw-r--r-- 1 user group 547 Apr 25 12:04 sim.R
[user@biowulf]$ **cat rfile1.R**
source("sim.R")
sim.fun(n.sim=50, output1="./rfile1.rds",seed=24963)
[user@biowulf]$ **cat rfile2.R**
source("sim.R")
sim.fun(n.sim=50, output1="./rfile2.rds",seed=27507)
[user@biowulf]$ **cat rfile.sw**
R --no-echo --no-restore --no-save < ./rfile1.R
R --no-echo --no-restore --no-save < ./rfile2.R
[user@biowulf]$ **swarm -f rfile.sw --time=10 --partition=quick --module R**
199110
[user@biowulf]$ **ls -lh \*.rds**
-rw-r--r-- 1 user group 445 Apr 25 12:52 rfile1.rds
-rw-r--r-- 1 user group 445 Apr 25 12:52 rfile2.rds

```

Full Rswarm usage:



```

Usage: Rswarm [options]
   --rfile=[file]   (required) R program requiring replication
   --sfile=[file]   (required) file with generated seeds, one per line
   --path=[path]    (required) directory for output of all files
   --reps=[i]       (required) number of replicates desired
   --sims=[i]       (required) number of sims per file
   --start=[i]      (required) starting file number
   --ext1=[string]    (optional) file extension for output file 1
   --ext2=[string]    (optional) file extension for output file 2`
   --help, -h         print this help text

```

Note that R scripts can be written to take a random seed as a command line argument
or derive it from the environment variable `SLURM_ARRAY_TASK_ID` to achieve
an equivalent result.



Using the parallel package
[top](#top)
 The R [parallel](http://stat.ethz.ch/R-manual/R-devel/library/parallel/doc/parallel.pdf) package provides functions for parallel execution of R code on
machines with multiple CPUs. Unlike other parallel processing methods, all jobs
share the full state of R when spawned, so no data or code needs to be
initialized if it was loaded before starting worker processes. The actual
spawning is very fast as well since no new R instance needs to be started.


#### Detecting the number of CPUs


Parallel includes the `dectectCores` function which is often used
to automatically detect the number of available CPUs. However, it always
reports *all* CPUs available on a node irrespective of how many CPUs
were allocated to the job. This is not the desired behavior for batch jobs
or sinteractive sessions. Instead, please use the `availableCores()` 
function from the future (or parallelly for R >= 4.0.3) package which correctly 
returns the number of *allocated* CPUs:



```

parallelly::availableCores() # for R >= 4.0.3
# or
future::availableCores()

```

Or, if you prefer, you could also write your own detection function
similar to the following example



```

detectBatchCPUs <- function() { 
    ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")) 
    if (is.na(ncores)) { 
        ncores <- as.integer(Sys.getenv("SLURM_JOB_CPUS_PER_NODE")) 
    } 
    if (is.na(ncores)) { 
        return(2)
    } 
    return(ncores) 
}

```

#### Random number generation


The state of the random number generator in each worker process has
to be carefully considered for any parallel workloads. See the help for
`mcparallel` and the parallel package documentation for
more details.


#### Example 1: mclapply



The `mclapply()` function calls `lapply()` in parallel, so that the first two arguments 
to `mclapply()` are exactly the same as for `lapply()`. Except the `mc.cores` argument needs 
to be specified to split the computatation across multiple CPUs on the same node. In most cases `mc.cores`
should be equal to the number of allocated CPUs.



```

> ncpus <- parallelly::availableCores()
> options(mc.cores = ncpus) # set a global option for parallel packages
# Then run mclapply() 
> mclapply(X, FUN, ..., mc.cores = ncpus)

```

Performance comparision between `lapply()` and `mclapply()`:



```

> library(parallel)
> ncpus <- parallelly::availableCores()
> N <- 10^6
> system.time(x<-lapply(1:N, function(i) {rnorm(300)}))
##   user  system elapsed
## 36.588   1.375  38.053
> system.time(x<-mclapply(1:N, function(i) {rnorm(300)},mc.cores = ncpus)) #Test on a phase5 node with ncpus=12
##   user  system elapsed
## 11.587  14.547  13.684

```

In this example, using 12 CPUs with `mclapply()` only reduced runtime by only 2.8 fold compared to running
on a single CPU. Under ideal conditions, the reduction would have been expected to be 12 fold. This means the work done
per CPU was less in the parallel case than in the sequential (single CPU) case. This is called *parallel efficiency*.
In this example the efficiency would have been sequential CPU time / parallel CPU time = (38.1 \* 1) / (13.7 \* 12) = 23%.


Parallel jobs should aim for an efficiency of 70-80%. Because parallel algorithms rarely
scale ideally to multiple CPUs we highly recommend performing scaling test before running programs in parallel.
To better optimize the usage of `mclapply()`, we benchmarked the performance of `mclappy()` with 2-32 CPUs 
and compared their efficiency:



![](/images/mclapply_benchmark.png)


The code used for benchmark was:



```

> library(parallel)
> library(microbenchmark)
> N <- 10^5
# benchmark the performance with 2-32 CPUs for 20 times
> for (n in c(2,4,6,8,16,32)) {
microbenchmark(mctest = mclapply(1:N, function(i) {rnorm(30000)},mc.cores=n),times = 20)
}

```

As show in the figure, this particular `mclapply()` should be run with no more than 6 CPUs
to ensure a higher than 70% of efficiency. This may be different for your code and should be tested
for each type of workload. Note that memory usage increases with more CPUs are used which
makes it even more important to not allocate more CPUs than necessary.


#### Example 2: foreach


A very convenient way to do parallel computations is provided by the 
[foreach](https://cran.r-project.org/web/packages/foreach/index.html) package. Here
is a simple example (copied from this 
[blog post](https://www.r-bloggers.com/the-wonders-of-foreach/))



```

> library(foreach)
> library(doParallel)
> library(doMC)
> registerDoMC(cores=future::availableCores())
> max.eig <- function(N, sigma) {
     d <- matrix(rnorm(N**2, sd = sigma), nrow = N)
     E <- eigen(d)$values
     abs(E)[[1]] } 
> library(rbenchmark)
> benchmark(
     foreach(n = 1:100) %do% max.eig(n, 1),
     foreach(n = 1:100) %dopar% max.eig(n, 1) )

##          test                             replications elapsed relative user.self sys.self user.child sys.child
##1    foreach(n = 1:100) %do% max.eig(n, 1)          100  32.696    3.243 32.632   0.059      0.000      0.00
##2 foreach(n = 1:100) %dopar% max.eig(n, 1)          100  10.083    1.000 3.037    3.389     43.417     10.73
>                                                                       

```

Note that with 12 CPUs we got a speedup of only 3.2 relative to sequential
resulting in a low parallel efficiency. Another cautionary tale to
carefully test scaling of parallel code.


A second way to run foreach in parallel:



```

> library(doParallel)
> cl <- makeCluster(future::availableCores())
> registerDoParallel(cl) 
 # parallel command
> ...
 # stop cluster
> stopCluster(cl)

```

What if we increased the number of tasks and the size of the largest matrix
(i.e. more work per task)? In the example above that means increasing the
`i` in `foreach(n=1:i)` using a fixed number of CPUs (32
in this case). We then calculated the speedup relative to execution on 2 CPUs:



.aligncenter {
 text-align: center;
}


![](/images/foreach_i_dopar-1.png)



If parallelism was 100% efficient, the speedup would be 16-fold. We recommend
running jobs at 70% parallel efficiency or better which would correspond to a
11-fold speedup in this case (blue horizontal line). In this example, 70% efficiency
is reached at `i > 300`. That means on biowulf you should only
run this code on 32 CPUs if for `i > 300`.


How does the code perform with different numbers of CPUs for `i = 500`.
Based on the results shown below, this code should be run with no more than 32 CPUs
to ensure that efficiency is better than 70%.



![](/images/foreach_i500_cpus_constrain-1.png)




Using the BiocParallel package
[top](#top)
 The R [BiocParallel](http://bioconductor.org/packages/release/bioc/html/BiocParallel.html) provides modified versions and novel implementation of
functions for parallel evaluation, tailored to use with Bioconductor objects. Like
the parallel package, it is not aware of slurm allocations and will therefore,
by default, try to use `parallel::detectCores() - 2` CPUs, which is
all but 2 CPUs installed on a compute node irrespective of how many CPUs have
been allocated to a job. That will lead to overloaded jobs and very inefficient
code. You can verify this by checking on the registered backends after allocating
an interactive session with 2 CPUs:



```

> library(BiocParallel)
> registered()
$MulticoreParam
class: MulticoreParam
  bpisup: FALSE; bpnworkers: 54; bptasks: 0; bpjobname: BPJOB
  bplog: FALSE; bpthreshold: INFO; bpstopOnError: TRUE
  bptimeout: 2592000; bpprogressbar: FALSE
  bpRNGseed: 
  bplogdir: NA
  bpresultdir: NA
  cluster type: FORK
[...snip...]

```

So the default backend (top of the registered stack) would use 54 workers on
2 CPUs. The default backend can be changed with



```

> options(MulticoreParam=quote(MulticoreParam(workers=future::availableCores())))
> registered()
$MulticoreParam
class: MulticoreParam
  bpisup: FALSE; bpnworkers: 2; bptasks: 0; bpjobname: BPJOB
  bplog: FALSE; bpthreshold: INFO; bpstopOnError: TRUE
[...snip..]

```

or



```

> register(MulticoreParam(workers = future::availableCores()), default=TRUE)

```

Alternatively, a param object can be passed to BiocParallel functions.



Implicit multithreading
[top](#top)
R can do implicit multithreading when using a subset of optimized functions
in the library or functions that take advantage of parallelized routines in the
lower level math libraries.
The function `crossprod(m)` which is equivalent to calculating
`t(m) %*% m`, for example, makes use of implicit parallelism in the
underlying math libraries and can benefit from using more than one thread. The number
of threads used by such functions is regulated by the environment variable
`OMP_NUM_THREADS`, which the R module sets automatically when
loaded as part of a batch or interactive job. Here is the runtime of
this function with different values for `OMP_NUM_THREADS`:


![crossprod benchmark](/images/R_implicitThreads2-mkl.png)

The code used for this benchmark was



```

# this file is benchmark2.R
runs <- 3
o <- 2^13
b <- 0

for (i in 1:runs) {
  a <- matrix(rnorm(o*o), o, o)
  invisible(gc())
  timing <- system.time({
    b <- crossprod(a)		# equivalent to: b <- t(a) %*% a
  })[3]
  cat(sprintf("%f\n", timing))
}

```

And was called with



```

node$ **module load R/3.5**
node$ **OMP\_NUM\_THREADS=1 Rscript benchmark2.R**
node$ **OMP\_NUM\_THREADS=2 Rscript benchmark2.R**
...
node$ **OMP\_NUM\_THREADS=32 Rscript benchmark2.R**

```

From within a job that had been allocated 32 CPUs.


Notes:


* The R module sets `OMP_NUM_THREADS` to 1.
 This is done to avoid problems when also using explicit parallelism that
 uses forked worker processes like for example the 'parallel' package. In that
 case each worker may start `OMP_NUM_THREADS` threads so there
 could be as many as `ncpus * OMP_NUM_THREADS` total threads. If 
 ncpus is set to the number of allocated CPUs, as would be the most common case,
 that would lead to a huge overload of the job and could potentially fail
 if the ulimit on the number of processes is exceeded.
* Increasing `OMP_NUM_THREADS` and allocating more CPUs
 may or may not increase performance of your code. Before running jobs
 with more CPUs it is vital to benchmark single jobs and prove that
 your code can benefit from implicit parallelism. Otherwise
 resources will be wasted.
* Benchmarking is also important to measure efficiency of parallelism. Clearly,
 in this example, running with more than 8 threads would be highly inefficient.


There appears to also be another level of parallelism within the R libraries.
One function that takes advantage of this is the `dist` function.
The level of parallelism allowed with this mechanism seems to be set with
two internal R functions (`setMaxNumMathThreads` and `setNumMathThreads`).
Note that this is a distinct mechanism - i.e. setting OMP\_NUM\_THREADS has
no impact on `dist` and `setMaxNumMathThreads`
has no impact on the performance of `crossprod`. Here is the performance
of `dist` with different numbers of threads:




![dist benchmark](/images/R_implicitThreads1-mkl.png)

The timings for this example were created with



```

# this file is benchmark1.R
rt <- data.frame()
o <- 2^12
m <- matrix(rnorm(o*o), o, o)
for (nt in c(1, 2, 4, 8, 16, 32)) {
    .Internal(setMaxNumMathThreads(nt)) 
    .Internal(setNumMathThreads(nt))
    res <- system.time(d <- dist(m))
    rt <- rbind(rt, c(nt, o, res[3]))
}
colnames(rt) <- c("threads", "order", "elapsed")
write.csv(rt, file="benchmark1.csv", row.names=F)

```

This was run within an allocation with 32 CPUs with



```

node$ **OMP\_NUM\_THREADS=1 Rscript benchmark1.R**

```

The same notes about benchmarking as above apply. Also note that there is
very little documentation about this to be found online.



R MPI jobs
[top](#top)
Our R installations include the [Rmpi](http://www.stats.uwo.ca/faculty/yu/Rmpi/) and
[pbdMPI](https://github.com/RBigData/pbdMPI) interfaces to MPI (OpenMPI in our case). R/MPI code
can be run as batch jobs or from an sinteractive session with `mpiexec` or `srun --mpi=pmix`.
Running MPI code from an interactive R session is currently not supported.


The higher level snow MPI cluster interface is currently not supported. However, 
the [doMPI](https://cran.r-project.org/web/packages/doMPI/vignettes/doMPI.pdf)
parallel backend for [foreach](https://cran.r-project.org/web/packages/foreach/index.html) is supported.


See our [MPI docs](/development/MPI.html) for more detail


### Example Rmpi code


This is a lower level Rmpi script



```

# this script is test1.r
library(Rmpi)
id <- mpi.comm.rank(comm=0)
np <- mpi.comm.size (comm=0)
hostname <- mpi.get.processor.name()

msg <- sprintf ("Hello world from task %03d of %03d, on host %s \n", id , np , hostname)
cat(msg)

invisible(mpi.barrier(comm=0))
invisible(mpi.finalize())

```


It can be submitted as a batch job with the following script:



```

#! /bin/bash
# this script is test1.sh

module load R/4.1.0 || exit 1
srun --mpi=pmix Rscript test1.r
## or
# mpiexec Rscript test1.r

```

which would be submitted with



```

[user@biowulf]$ **sbatch --ntasks=4 --nodes=2 --partition=multinode test1.sh**

```

And would generate output similar to



```

Hello world from task 000 of 004, on host cn4101
Hello world from task 001 of 004, on host cn4102
Hello world from task 002 of 004, on host cn4103
Hello world from task 003 of 004, on host cn4104

```

Here is a Rmpi example with actual collective communication though still very simplistic.
This script derives an estimate for π in each task, gathers the results in task 0 and
repeats this process n times to arrive at a final estimate:



```

# this is test2.r
library(Rmpi)

# return a random number from /dev/urandom
readRandom <- function() {
  dev <- "/dev/urandom"
  rng = file(dev,"rb", raw=TRUE)
  n = readBin(rng, what="integer", 1) # read some 8-byte integers 
  close(rng)
  return( n[1] ) # reduce range and shift
}

pi_dart <- function(i) {
    est <- mean(sqrt(runif(throws)^2 + runif(throws)^2) <= 1) * 4
    return(est)
}

id <- mpi.comm.rank(comm=0)
np <- mpi.comm.size (comm=0)
hostname <- mpi.get.processor.name()
rngseed <- readRandom()
cat(sprintf("This is task %03d of %03d, on host %s with seed %i\n", 
    id , np , hostname, rngseed))
set.seed(rngseed)

throws <- 1e7
rounds <- 400
pi_global_sum = 0.0
for (i in 1:rounds) {
    # each task comes up with its own estimate of pi
    pi_est <- mean(sqrt(runif(throws)^2 + runif(throws)^2) <= 1) * 4
    # then we gather them all in task 0; type=2 means that the values are doubles
    pi_task_sum <- mpi.reduce(pi_est, type=2, op="sum", dest=0, comm=0)
    if (id == 0) {
        # task with id 0 then uses the sum to calculate an avarage across the
        # tasks and adds that to the global sum
        pi_global_sum <- pi_global_sum + (pi_task_sum / np)
    }
}

# when we're done, the task id 0 averages across all the rounds and prints the result
if (id == 0) {
    cat(sprintf("Real value of pi = %.10f\n", pi))
    cat(sprintf("  Estimate of pi = %.10f\n", pi_global_sum / rounds))
}

invisible(mpi.finalize())

```

Submitting this script with a batch script similar to the first example results
in output like this:



```

This is task 000 of 004, on host cn4101 with seed -303950071
This is task 001 of 004, on host cn4102 with seed -1074523673
This is task 002 of 004, on host cn4103 with seed 788983269
This is task 003 of 004, on host cn4104 with seed -922785662
Real value of pi = 3.1415926536
  Estimate of pi = 3.1415935438

```

### doMPI example code


doMPI provides an MPI backend for the foreach package. Here is a simple hello-world-ish
doMPI example. Note that in my testing the least issues were encountered when
the foreach loops were run from the first rank of the MPI job.



```

suppressMessages({
    library(Rmpi)
    library(doMPI)
    library(foreach)
})

myrank <- mpi.comm.rank(comm=0)
cl <- doMPI::startMPIcluster()
registerDoMPI(cl)
if (myrank == 0) {
    cat("-------------------------------------------------\n")
    cat("== This is rank 0 running the foreach loops ==\n")
    

    x <- foreach(i=1:16, .combine="c") %dopar% {
        id <- mpi.comm.rank(comm=0)
        np <- mpi.comm.size (comm=0)
        hostname <- mpi.get.processor.name()
        sprintf ("Hello world from process %03d of %03d, on host %s \n", id , np , hostname)
    }
    print(x)
    
    x <- foreach(i=1:200, .combine="c") %dopar% {
        sqrt(i)
    }

    cat("-------------------------------------------------\n")
    print(x)

    x <- foreach(i=1:16, .combine="cbind") %dopar% {
        set.seed(i)
        rnorm(3)
    }


    cat("-------------------------------------------------\n")
    print(x)

}
closeCluster(cl)
mpi.quit(save="no")

```

Do not specify the number of tasks to run. doMPI clusters will hang during shutdown
if doMPI has to spawn worker processes. Instead, let `mpiexec` start the processes and then
doMPI will wire up a main process and workers from the existing processes. Note that the 
`startMPIcluster` function has to be called early in the script for that reason.



Run a shiny app on biowulf
[top](#top)
Shiny apps can be run on biowulf for a single user. Since they require tunneling a running
shiny app cannot be shared with other users. However, if the code for the app is accessible
to different users they each can run an ephemeral shiny app on biowulf. We will use the following
example application:



```

## this file is 01_hello.r
library(shiny)


ui <- fluidPage(

  # App title ----
  titlePanel("Hello Shiny!"),
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    # Sidebar panel for inputs ----
    sidebarPanel(
      # Input: Slider for the number of bins ----
      sliderInput(inputId = "bins",
                  label = "Number of bins:",
                  min = 1,
                  max = 50,
                  value = 30)
    ),
    # Main panel for displaying outputs ----
    mainPanel(
      # Output: Histogram ----
      plotOutput(outputId = "distPlot")
    )
  )
)

# Define server logic required to draw a histogram ----
server <- function(input, output) {
  # Histogram of the Old Faithful Geyser Data ----
  # with requested number of bins
  # This expression that generates a histogram is wrapped in a call
  # to renderPlot to indicate that:
  #
  # 1. It is "reactive" and therefore should be automatically
  #    re-executed when inputs (input$bins) change
  # 2. Its output type is a plot
  output$distPlot <- renderPlot({
    x    <- faithful$waiting
    bins <- seq(min(x), max(x), length.out = input$bins + 1)
    hist(x, breaks = bins, col = "#007bc2", border = "white",
         xlab = "Waiting time to next eruption (in mins)",
         main = "Histogram of waiting times")

    })
}

# which port to run on
port <- tryCatch(
  as.integer(Sys.getenv("PORT1", "none")),
  error = function(e) {
    cat("Please remember to use --tunnel to run a shiny app")
    cat("See https://hpc.nih.gov/docs/tunneling/")
    stop()
  }
)

# run the app
shinyApp(
  ui,
  server,
  options = list(port=port, launch.browser=F, host="127.0.0.1")
)

```

Start an sinteractive session with a [tunnel](https://hpc.nih.gov/docs/tunneling/)



```

[user@biowulf]$ **sinteractive --cpus-per-task=2 --mem=6g --gres=lscratch:10 --tunnel**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job
[user@cn3144 ~]$ **module load R**
[user@cn3144 ~]$ **Rscript 01\_hello.r**

Listening on http://127.0.0.1:34239

```

After you set up your tunnel you can use the URL above to access the shiny app.



Notes for individual packages
[top](#top)
### h2o


h2o is a machine learning package written in java. The R interface starts a java h2o instance
with a given number of threads and then connects to it through http. This fails on compute nodes
if the http proxy variables are set. Therefore it is necessary to unset `http_proxy`
before using h2o:



```

[user@biowulf]$ **sinteractive --cpus-per-task=4 --mem=20g --gres=lscratch:10**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job
[user@cn3144 ~]$ **module load R/4.2**
[user@cn3144 ~]$ **unset http\_proxy**
[user@cn3144 ~]$ **R**
R version 4.2.0 (2022-04-22) -- "Vigorous Calisthenics"
Copyright (C) 2022 The R Foundation for Statistical Computing
Platform: x86_64-pc-linux-gnu (64-bit)

[...snip...]
> library(h2o)
> h2o.init(ip='localhost', nthreads=future::availableCores(), max_mem_size='12g')
H2O is not running yet, starting it now...

Note:  In case of errors look at the following log files:
    /tmp/RtmpVdW92Y/h2o_user_started_from_r.out
    /tmp/RtmpVdW92Y/h2o_user_started_from_r.err

openjdk version "1.8.0_161"
OpenJDK Runtime Environment (build 1.8.0_161-b14)
OpenJDK 64-Bit Server VM (build 25.161-b14, mixed mode)

Starting H2O JVM and connecting: . Connection successful!

R is connected to the H2O cluster:
    H2O cluster uptime:         1 seconds 683 milliseconds
    H2O cluster timezone:       America/New_York
    H2O data parsing timezone:  UTC
    H2O cluster version:        3.36.1.2
    H2O cluster version age:    3 months and 20 days !!!
    H2O cluster name:           H2O_started_from_R_user_ywu882
    H2O cluster total nodes:    1
    H2O cluster total memory:   10.64 GB
    H2O cluster total cores:    4
    H2O cluster allowed cores:  4
    H2O cluster healthy:        TRUE
    H2O Connection ip:          localhost
    H2O Connection port:        54321
    H2O Connection proxy:       NA
    H2O Internal Security:      FALSE
    R Version:                  R version 4.2.0 (2022-04-22)

>

```

### dyno


dyno is a meta package that installs several other packages from the dynvers
(<https://github.com/dynverse>). It
includes some cran packages and some packages only available on github. We
generally don't install any new github-only R packages any more so here are the
instructions for installing this as a user.


##### Installation



```

###
### 1. install with the default dependent packages
###
[user@biowulf]$ **sinteractive --gres=lscratch:5**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job
[user@cn3144 ~]$ **module load R/4.2**
[user@cn3144 ~]$ **R -q --no-save --no-restore -e 'devtools::install\_github("dynverse/dyno")'**

###
### 2. install a pached version of babelwhale
###
[user@cn3144 ~]$ **git clone https://github.com/dynverse/babelwhale.git**
[user@cn3144 ~]$ **patch -p0 <<'\_\_EOF\_\_'
--- babelwhale/R/run.R.orig 2021-07-16 20:58:26.563714000 -0400
+++ babelwhale/R/run.R 2021-07-16 20:58:26.483721000 -0400
@@ -122,6 +122,8 @@
 environment\_variables %>% gsub("^.\*=", "", .),
 environment\_variables %>% gsub("^(.\*)=.\*$", "SINGULARITYENV\_\\1", .)
 ),
+ "http\_proxy" = Sys.getenv("http\_proxy"),
+ "https\_proxy" = Sys.getenv("https\_proxy"),
 "SINGULARITY\_TMPDIR" = tmpdir,
 "SINGULARITY\_CACHEDIR" = config$cache\_dir,
 "PATH" = Sys.getenv("PATH") # pass the path along
\_\_EOF\_\_**

[user@cn3144 ~]$ **R CMD INSTALL babelwhale**

###
### 3. Create a configuration that uses a dedicated singularity cache somewhere
###    outside the home directory. In this example using `/data/$USER/dynocache`
###

[user@cn3144 ~]$ **R -q --no-save --no-restore <<'\_\_EOF\_\_'
config <- babelwhale::create\_singularity\_config(
 cache\_dir = "/data/$USER/dynocache"
)
babelwhale::set\_default\_config(config, permanent = TRUE)
\_\_EOF\_\_**

```

##### Notes:


* There is a bug in the `plot\_dimred` function - see
 https://github.com/dynverse/dynplot/issues/54. I tried installing
 the devel version but it didn't help. YMMV
* If you want to use the parts that use shiny (e.g. guidelines\_shiny)
 you will need to use a tunnel and call the function like
 `guidelines_shiny(dataset, port=####, launch.browser=F)` where
 #### is the port set up by the `sinteractive --tunnel`



Documentation
[top](#top)
* [The R Homepage](http://www.cran.r-project.org/)
* [PDF manuals](http://cran.r-project.org/manuals.html) including an Introduction to R.
* [The R FAQ](http://cran.r-project.org/doc/FAQ/R-FAQ.html)
* ['parallel' documentation](http://stat.ethz.ch/R-manual/R-devel/library/parallel/doc/parallel.pdf)
* [Rmpi documentation](http://www.stats.uwo.ca/faculty/yu/Rmpi/)
* [Snow user guide](http://www.sfu.ca/~sblay/R/snow.html)


 
$(document).ready(function() { $('#r-packages').DataTable({
 columnDefs: [
 // Left align the header content of column 1
 { className: "dt-left", targets: [ '\_all' ] }
 ]
});});









