

document.querySelector('title').textContent = 'xtail: genome-wide assessment of differential translations with ribosome profiling data';
**xtail: genome-wide assessment of differential translations with ribosome profiling data
.**


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



Xtail is an analysis pipeline tailored for ribosome profiling data 
that comprehensively and accurately identifies differentially translated genes in pairwise comparisons. 
Applied on simulated and real datasets, Xtail exhibits high sensitivity with minimal false-positive rates, 
outperforming existing methods in the accuracy of quantifying differential translations.



### References:


* Zhengtao Xiao, Qin Zou, Yu Liu & Xuerui Yang   

*Genome-wide assessment of differential translations with ribosome profiling data*    

[Nature Communications volume 7, Article number: 11194 (2016)](https://www.nature.com/articles/ncomms11194)


Documentation
* [xtail GitHub page](https://github.com/xryanglab/xtail)
* [xtail Documentation](https://rdrr.io/github/xryanglab/xtail/man/)


Important Notes
* Module Name: xtail (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set
	+ **XTAIL\_HOME**  installation directory
	+ **XTAIL\_BIN**       executable directory
	+ **XTAIL\_SRC**       source code directory
	+ **XTAIL\_DATA**  sample data directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=4g**
[user@cn0911 ~]$**module load xtail** 
[+] Loading singularity  3.10.5  on cn4185
[+] Loading xtail  1.1.5
[user@cn0911 ~]$**wget https://github.com/xryanglab/xtail/archive/refs/tags/v1.1.4.tar.gz**
[user@cn0911 ~]$**tar -zxf v1.1.4.tar.gz && rm -f v1.1.4.tar.gz && cd xtail-1.1.4**
[user@cn0911 ~]$**R**

R version 3.6.1 (2019-07-05) -- "Action of the Toes"
...
> **library(xtail)** 
...
> **data(xtaildata)**

```

Get the mrna count data and rpf count data. For the example only the first
# 100 are used:

```

> **test.mrna <- xtaildata$mrna[1:100,]**
> **test.mrna** 
                control1 control2 treat1 treat2
ENSG00000000003      825      955    866   1039
ENSG00000000419     1054      967    992    888
ENSG00000000457       71       75    139     95
ENSG00000000460      191      162    199    201
ENSG00000000971       81        2     88     11
...
ENSG00000006530      562      588    547    597
ENSG00000006534      176      106    157    141
ENSG00000006555       32       36     40     52
> **test.rpf <- xtaildata$rpf[1:100,]**
> **test.rpf**
                control1 control2 treat1 treat2
ENSG00000000003      143      302    197    195
ENSG00000000419      234      481    383    306
ENSG00000000457       12       17     17     15
...
ENSG00000006530      127      214    130    176
ENSG00000006534       13       25     17     10
ENSG00000006555        3       11      4      3

```

Assign condition labels to samples:

```

> **condition <- c("control","control","treat","treat")**

```

Run xtail

```

> **test.results <- xtail(test.mrna,test.rpf,condition, threads = 2)**
Calculating the library size factors
1. Estimate the log2 fold change in mrna
2. Estimate the log2 fold change in rpf
3. Estimate the difference between two log2 fold changes
4. Estimate the log2 ratio in first condition
5. Estimate the log2 ratio in second condition
6. Estimate the difference between two log2 ratios
Number of the log2FC and log2R used in determining the final p-value
 log2FC: 16
 log2R: 84
>  **q()**
[user@cn0911 ~]$**exit**

```
 
End the interactive session:

```

[user@cn0911 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```





