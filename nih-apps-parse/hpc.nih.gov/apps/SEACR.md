

document.querySelector('title').textContent = 'SEACR: Sparse Enrichment Analysis for CUT&RUN';
**SEACR: Sparse Enrichment Analysis for CUT&RUN**


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



SEACR is intended to call peaks and enriched regions 
from sparse Cleavage Under Targets and Release Using Nuclease (CUT&RUN)
or chromatin profiling data 
in which background is dominated by "zeroes" (i.e. regions with no read coverage). 



### References:


* Michael P. Meers, Terri Bryson, Steven Henikoff,   

"A streamlined protocol and analysis pipeline for CUT&RUN chromatin profiling",
  

bioRxiv, 2019, doi: https://doi.org/10.1101/569129.
* Peter J. Skene and Steven Henikoff,   

"An efficient targeted nuclease strategy for high-resolution mapping of DNA binding sites",   

eLife, 2017, Jan 16. doi: 10.7554/eLife.21856.


Documentation
	+ [SEACR GitHub Page](https://github.com/mpmeers/SEACR)Important Notes
	+ Module Name: SEACR (see [the modules page](/apps/modules.html) for more information)
	+ Unusual environment variables set
		- **SEACR\_HOME** installation directory
		- **SEACR\_BIN** executable folder
		- **SEACR\_DATA** sample data for running SEACR
Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=16g**
[user@cn3144 ~]$ **module load SEACR** 
[+] Loading GSL 2.6 for GCC 9.2.0 ...
[+] Loading gcc  9.2.0  ...
[+] Loading openmpi 3.1.4  for GCC 9.2.0
[+] Loading ImageMagick  7.0.8  on cn2357
[+] Loading HDF5  1.10.4
[+] Loading pandoc  2.9.1  on cn2357
[+] Loading R 3.6.1
[+] Loading bedtools  2.29.0
[+] Loading SEACR  1.3

```

Display the usage message for SEACR:

```

[user@cn3144 ~]$ **SEACR.sh** 
	SEACR: Sparse Enrichment Analysis for CUT&RUN
	
	Usage: bash SEACR.sh .bg [.bg | ] [norm | non] [union | AUC] output prefix
 
 Description of input fields:
 
 Field 1: Target data bedgraph file in UCSC bedgraph format (https://genome.ucsc.edu/goldenpath/help/bedgraph.html) that omits regions containing 0 signal.
 
 Field 2: Control (IgG) data bedgraph file to generate an empirical threshold for peak calling. Alternatively, a numeric threshold n between 0 and 1 returns the top n fraction of peaks based on total signal within peaks.
 
 Field 3: “norm” denotes normalization of control to target data, “non” skips this behavior. norm is recommended unless experimental and control data are already rigorously normalized to each other (e.g. via spike-in).
 
 Field 4: “union” forces implementation of a maximum signal threshold in addition to the total signal threshold, and corresponds to the “union” mode described in the text, whereas “AUC” avoids this behavior, and corresponds to “AUC only” mode.
 
 Field 5: Output prefix
 
 Output file:

 .auc.threshold.merge.bed (Bed file of enriched regions)
 
 Output data structure: 
 
      
 
 Description of output fields:

 Field 1: Chromosome
 
 Field 2: Start coordinate
 
 Field 3: End coordinate
 
 Field 4: Total signal contained within denoted coordinates
 
 Field 5: Maximum bedgraph signal attained at any base pair within denoted coordinates
 
 Field 6: Region representing the farthest upstream and farthest downstream bases within the denoted coordinates that are represented by the maximum bedgraph signal
 
 Examples:

 bash SEACR.sh target.bedgraph IgG.bedgraph norm AUC output
 Calls enriched regions in target data using normalized IgG control track with AUC threshold
 
 bash SEACR.sh target.bedgraph IgG.bedgraph non union output
 Calls enriched regions in target data using non-normalized IgG control track with AUC and max signal thresholds

 bash SEACR.sh target.bedgraph 0.01 non AUC output
 Calls enriched regions in target data by selecting the top 1% of regions by area under the curve (AUC)


```

Copy sample data from an application folder to your current folder:

```

[user@cn3144 ~]$ **cp $SEACR\_DATA/\* .** 

```

This command will copy the following six files: 

```

[user@cn3144 ~]$ **ls**
DE_FoxA2.hg19.bed  DE_IgG.hg19.bed  DE_Sox2.hg19.bed  
H1_FoxA2.hg19.bed  H1_IgG.hg19.bed  H1_Sox2.hg19.bed

```

The data have been downloaded from the [GEO website](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126612). Here,   

- IgG stands for immunoglobulin G (IgG); IgG experiments are used as control   

- Sox2 and FoxA2 are two transcription factors in   

 - human embryonic stem cells (hESCs) and  

 - Definitive Endoderm (DE) cells  

- Sox2 expression is restricted to hESCs cells, and  

- FoxA2 expression is restricted to DE cells.   

  

Prepare files in the bedGraph format, which will be used as inputs for SEACR:

```

[user@cn3144 ~]$ **ln -s /fdb/igenomes/Homo\_sapiens/UCSC/hg19/GenomeStudio/Homo\_sapiens/UCSC-hg19/ChromInfo.txt**
[user@cn3144 ~]$ **genomeCoverageBed -i DE\_FoxA2.hg19.bed -g ChromInfo.txt -bg > DE\_FoxA2.hg19.bedGraph**
[user@cn3144 ~]$ **genomeCoverageBed -i H1\_FoxA2.hg19.bed -g ChromInfo.txt -bg > H1\_FoxA2.hg19.bedGraph**

[user@cn3144 ~]$ **genomeCoverageBed -i DE\_IgG.hg19.bed -g ChromInfo.txt -bg > DE\_IgG.hg19.bedGraph**
[user@cn3144 ~]$ **genomeCoverageBed -i H1\_IgG.hg19.bed -g ChromInfo.txt -bg > H1\_IgG.hg19.bedGraph**

[user@cn3144 ~]$ **genomeCoverageBed -i DE\_Sox2.hg19.bed -g ChromInfo.txt -bg > DE\_SoxA2.hg19.bedGraph**
[user@cn3144 ~]$ **genomeCoverageBed -i H1\_Sox2.hg19.bed -g ChromInfo.txt -bg > H1\_SoxA2.hg19.bedGraph**

```

Call enriched regions in target data using normalized IgG control track with AUC threshold:

```

[user@cn3144 ~]$ **SEACR.sh DE\_FoxA2.hg19.bedGraph DE\_IgG.hg19.bedGraph norm stringent AUC output**
Calling enriched regions with control file
Must specify "norm" for normalized or "non" for non-normalized data processing in third input
[helixapp@cn4470 SEACR]$  SEACR.sh DE_FoxA2.hg19.bedGraph DE_IgG.hg19.bedGraph norm stringent AUC output
Calling enriched regions with control file
Normalizing control to experimental bedgraph
Using stringent threshold
Creating experimental AUC file: Wed Jan 29 13:55:19 EST 2020
Creating control AUC file: Wed Jan 29 13:55:49 EST 2020
Calculating optimal AUC threshold: Wed Jan 29 13:56:17 EST 2020
Calculating threshold using normalized control: Wed Jan 29 13:56:17 EST 2020
Warning messages:
1: In ((pctremain(x0) + min(pctremain(z)))/2) - pctremain(z) :
  longer object length is not a multiple of shorter object length
2: In ((pctremain(x0) + min(pctremain(z)))/2) - pctremain(z) :
  longer object length is not a multiple of shorter object length
Creating thresholded feature file: Wed Jan 29 13:58:05 EST 2020
Empirical false discovery rate = 0.356839643073724
Merging nearby features and eliminating control-enriched features: Wed Jan 29 13:58:20 EST 2020
Removing temporary files: Wed Jan 29 13:58:20 EST 2020
Done: Wed Jan 29 13:58:20 EST 2020

```

Likewise,

```

[user@cn3144 ~]$ **SEACR.sh DE\_Sox2.hg19.bedGraph DE\_IgG.hg19.bedGraph norm stringent AUC output\_DE\_Sox2\_norm\_IgG**
...
[user@cn3144 ~]$ **SEACR.sh H1\_FoxA2.hg19.bedGraph H1\_IgG.hg19.bedGraph norm stringent AUC output\_H1\_FoxA2\_norm\_IgG**
...
[user@cn3144 ~]$ **SEACR.sh H1\_Sox2.hg19.bedGraph H1\_IgG.hg19.bedGraph norm stringent AUC output\_H1\_Sox2\_norm\_IgG**
...

```

Now, call enriched regions in target data in a different mode, using non-normalized IgG control track with AUC and max signal thresholds:

```

[user@cn3144 ~]$ **SEACR.sh DE\_FoxA2.hg19.bedGraph DE\_IgG.hg19.bedGraph non stringent union output\_DE\_FoxA2\_non-norm\_IgG**
...
[user@cn3144 ~]$ **SEACR.sh DE\_Sox2.hg19.bedGraph DE\_IgG.hg19.bedGraph non stringent union output\_DE\_Sox2\_non-norm\_IgG**
...
[user@cn3144 ~]$ **SEACR.sh H1\_FoxA2.hg19.bedGraph H1\_IgG.hg19.bedGraph non stringent union output\_H1\_FoxA2\_non-norm\_IgG**
...
[user@cn3144 ~]$ **SEACR.sh H1\_Sox2.hg19.bedGraph H1\_IgG.hg19.bedGraph non stringent union output\_H1\_Sox2\_non-norm\_IgG**
...

```

As expected, the output indicates that in both the modes of runn ing the software, SEACR detected a large number of called peaks in the cell types where the Sox2 factor or FoxA2 factor is expressed:   


```

H1_Sox2_norm_IgG      - 17047 peaks
H1_Sox2_non-norm_IgG  - 20267 peaks 
DE_FoxA2_norm_IgG     - 6227  peaks 
DE_FoxA2_non-norm_IgG - 7436 peaks 

```

and a small number of (spurious) peaks in the cells where the corresponding factor is not expressed:   


```

H1_FoxA2_norm_IgG     - 1 peak 
H1_FoxA2_non-norm_IgG - 4 peaks 
DE_Sox2_norm_IgG      - 2 peaks 
DE_Sox2_non-norm_IgG  - 4 peaks 

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. SEACR.sh). For example:



```

#!/bin/bash
module load SEACR     

genomeCoverageBed -i DE_IgG.hg19.bed -g ChromInfo.txt -bg   > DE_IgG.hg19.bedGraph
genomeCoverageBed -i H1_IgG.hg19.bed -g ChromInfo.txt -bg   > H1_IgG.hg19.bedGraph
genomeCoverageBed -i DE_FoxA2.hg19.bed -g ChromInfo.txt -bg > DE_FoxA2.hg19.bedGraph
genomeCoverageBed -i H1_FoxA2.hg19.bed -g ChromInfo.txt -bg > H1_FoxA2.hg19.bedGraph
genomeCoverageBed -i DE_Sox2.hg19.bed -g ChromInfo.txt -bg  > DE_Sox2.hg19.bedGraph
genomeCoverageBed -i H1_Sox2.hg19.bed -g ChromInfo.txt -bg  > H1_Sox2.hg19.bedGraph

SEACR.sh DE_FoxA2.hg19.bedGraph DE_IgG.hg19.bedGraph norm AUC  output_DE_FoxA2_norm_IgG
SEACR.sh DE_Sox2.hg19.bedGraph  DE_IgG.hg19.bedGraph norm AUC  output_DE_Sox2_norm_IgG
SEACR.sh H1_FoxA2.hg19.bedGraph H1_IgG.hg19.bedGraph norm AUC  output_H1_FoxA2_norm_IgG
SEACR.sh H1_Sox2.hg19.bedGraph  H1_IgG.hg19.bedGraph norm AUC  output_H1_Sox2_norm_IgG
...
etc.


```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] SEACR.sh
```
