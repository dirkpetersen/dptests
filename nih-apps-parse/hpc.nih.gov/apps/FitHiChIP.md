

document.querySelector('title').textContent = 'FitHiChIP: Identification of significant chromatin contacts from HiChIP data';
**FitHiChIP: Identification of significant chromatin contacts from HiChIP data**


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



FitHiChIP is a computational method for identifying chromatin contacts 
among regulatory regions such as enhancers and promoters from HiChIP/PLAC-seq data. 
FitHiChIP jointly models the non-uniform coverage and genomic distance scaling of HiChIP data,
captures previously validated enhancer interactions for several genes 
including MYC and TP53, and recovers contacts genome-wide 
that are supported by ChIA-PET, promoter capture Hi-C and Hi-C data.



### References:


* Sourya Bhattacharyya, Vivek Chandra, Pandurangan Vijayanand, and Ferhat Ay,   

*Identification of significant chromatin contacts from HiChIP data by FitHiChIP.*   

[*Nature Communications* volume 10, Article number: 4221 (2019) .](https://www.nature.com/articles/s41467-019-11950-y)


Documentation
* [FitHiChIP Github page](https://github.com/ay-lab/FitHiChIP)
* [FitHiChIP Manual](https://ay-lab.github.io/FitHiChIP)


Important Notes
* Module Name: FitHiChIP (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Implemented as a Singularity container
* Unusual environment variables set
	+ **FITHICHIP\_HOME**  installation directory
	+ **FITHICHIP\_BIN**    executables directory
	+ **FITHICHIP\_SRC**    source directory
	+ **FITHICHIP\_DATA**  sample data directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=8g -c10 --gres=lscratch:20**
[user@cn3316 ~]$ **module load FitHiChIP** 
[+] Loading singularity/3.10.5  on cn3316
[+] Loading FitHiChIP  11.0
[user@cn3316 ~]$ **mkdir -p /data/$USER/FitHiChIP && cd /data/$USER/FitHiChIP**

```

Copy the source code folder and the sample sample data folder to your current directory: 

```

[user@cn3316 ~]$ **cp -r $FITHICHIP\_SRC/\* .**
[user@cn3316 ~]$ **cp -r $FITHICHIP\_DATA .**

```

Run the FitHiChIP\_HiCPro.sh executable in the source folder: 

```

[user@cn3316 ~]$ **./FitHiChIP\_HiCPro.sh -C configfile\_test**
 
 ================ Parsing input configuration file =================
...

 ================ Verifying input configuration parameters =================
...

 ====== Changing relative pathnames of the input files to their absolute path names ==========
...

 ====== Final checing of input parameters ==========

 ====== Writing input parameters ==========
...

 ================ Processing HiC-pro generated valid pairs and / or matrix files provided as input =================

 ====>> Computing HiC-pro matrices from the input valid pairs file

 ====>> Executable to generate contact matrix from valid pairs: /usr/local/apps/FitHiChIP/11.0/bin/build_matrix

***** HiC-pro input valid pairs file in gzipped format

 ====>> Created file : /data/user/FitHiChIP/TestData/results/Matrix_BinSize5000/Matrix_abs.bed

 ====>> Created file : /data/user/FitHiChIP/TestData/results/Matrix_BinSize5000/Matrix.matrix

 ================ Creating input interactions (bin pairs + CC) =================
...

 ======= Generated interaction file : /data/user/FitHiChIP/TestData/results/Matrix_BinSize5000/FitHiChIP.interactions.initial.bed
==>>> Number of locus pairs with nonzero contact count (without any distance thresholding): 2364391


 ======= Limiting input interactions to the specified distance ranges 5000 to 1000000 =========
...

===>> Number of cis pairs with nonzero contact count (after distance thresholding): 726467


 ================ Generating coverage statistics and bias for individual bins =================
...
 Computing 1D coverage - Processing chromosome : chr1  -- number of input peaks for this chromosome : 4175
 Computing 1D coverage - Processing chromosome : chr10  -- number of input peaks for this chromosome : 1797
 Computing 1D coverage - Processing chromosome : chr11  -- number of input peaks for this chromosome : 2155
...
 Computing 1D coverage - Processing chromosome : chrX  -- number of input peaks for this chromosome : 1161
 Computing 1D coverage - Processing chromosome : chrY  -- number of input peaks for this chromosome : 156
...
 ================ Computing bias statistics - coverage bias will be employed =================
R...

 ================ creating full feature file for FitHiChIP =================


======== Created full feature file : /data/user/FitHiChIP/TestData/results/NormFeatures/Coverage_Bias/FitHiChIP.AllBin_CompleteFeat.bed

 ================ Generating interactions + features for significance estimation =================
...

 **** Start of while Loop ----- current interaction type: 3  ******

 ============ Calling significant interactions ===============
...
===>> Total Number of input interactions (locus pairs): 312536
****** Total number of training interactions: 312536 ********
****** Number of contacts per bin (allowed for equal occupancy binning): 2746 ********
...
==>>> modeled the bias regression -- time: 0.646705150604248


 *** NumElem: 312536 ExpTotContact : 432569
 *** Modeling regression bias probability based P value --- uniq distance idx: 1  si: 1  ei: 26635 num elem : 26635
 *** Modeling regression bias probability based P value --- uniq distance idx: 2  si: 26636  ei: 45246 num elem : 18611
 *** Modeling regression bias probability based P value --- uniq distance idx: 3  si: 45247  ei: 59613 num elem : 14367
...
 *** Modeling regression bias probability based P value --- uniq distance idx: 198  si: 311681  ei: 311962 num elem : 282
 *** Modeling regression bias probability based P value --- uniq distance idx: 199  si: 311963  ei: 312238 num elem : 276
 *** Modeling regression bias probability based P value --- uniq distance idx: 200  si: 312239  ei: 312536 num elem : 298

 *** p-value estimation is complete for the bias regression ****

******** FINISHED calling significant interactions
----- Extracted significant interactions ---- FDR threshold lower than: 0.01
...
 ********** Merged filtering option is true ************


******** applying merge filtering on the FitHiChIP significant interactions ******

****** Merge filtering of adjacent loops is enabled *****
***** within function of merged filtering - printing the parameters ***
*** bin_size:  5000
*** headerInp:  1
*** connectivity_rule:  8
*** TopPctElem:  100
*** NeighborHoodBinThr:  10000
*** QValCol:  26
*** PValCol:  25
*** SortOrder:  0
==================== End of merge filtering adjacent interactions !!! ======================
----- Applied merged filtering (connected component model) on the adjacent loops of FitHiChIP
pass1 - making usageList (24 chroms): 1 millis
pass2 - checking and writing primary data (1882 records, 18 fields): 13 millis
/bin/bash: /opt/conda/envs/fithichip/lib/libtinfo.so.6: no version information available (required by /bin/bash)

 ---- Within function of plotting distance vs contact count ----
 Input interaction file: /data/user/FitHiChIP/TestData/results/FitHiChIP_Peak2ALL_b5000_L5000_U1000000/P2PBckgr_0/Coverage_Bias/FitHiC_BiasCorr/Merge_Nearby_Interactions/FitHiChIP.interactions_FitHiC_Q0.01_MergeNearContacts.bed
 Output plot file: /data/user/FitHiChIP/TestData/results/FitHiChIP_Peak2ALL_b5000_L5000_U1000000/P2PBckgr_0/Coverage_Bias/FitHiC_BiasCorr/Merge_Nearby_Interactions/FitHiChIP.interactions_FitHiC_Q0.01_MergeNearContacts_Dist_CC.png Warning message:
Removed 1 rows containing missing values (`geom_bar()`).
Merged filtering significant interactions - created various epigenome browser compatible files for these interactions!!!
Updated CurrIntType: 4

 **** Now summarizing FitHiChIP results in the HTML file ***

***** FitHiChIP pipeline is completely executed - congratulations !!! *****

```

Likewise, FitHiChIP can be run with other configuration files. 

```

[user@cn3316 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```

Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. fithichip.sh). For example:



```

#!/bin/bash
module load FitHiChIP
FitHiChIP_HiCPro.sh -C configfile_1
FitHiChIP_HiCPro.sh -C configfile_2
FitHiChIP_HiCPro.sh -C configfile_3
FitHiChIP_HiCPro.sh -C configfile_4

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] fithichip.sh
```





