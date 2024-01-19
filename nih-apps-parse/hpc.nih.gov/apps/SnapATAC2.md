

document.querySelector('title').textContent = 'SnapATAC2: Single Nucleus Analysis Pipeline for ATAC-seq';
**SnapATAC2: Single Nucleus Analysis Pipeline for ATAC-seq**


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



SnapATAC is a software package for analyzing scATAC-seq datasets. 
Conventional assays to map regulatory elements via open chromatin analysis of primary tissues 
is hindered by sample heterogeneity. 
SnapATAC dissects cellular heterogeneity in an unbiased manner and map the trajectories of cellular states.



### References:


* Rongxin Fang, Sebastian Preissl, Yang Li, Xiaomeng Hou, Jacinta Lucero, Xinxin Wang, Amir Motamedi, Andrew K. Shiau, Xinzhu Zhou, Fangming Xie, Eran A. Mukamel, Kai Zhang, Yanxiao Zhang, M. Margarita Behrens, Joseph R. Ecker & Bing Ren   

 *Comprehensive analysis of single cell ATAC-seq data with SnapATAC*   

[NATURE COMMUNICATIONS | (2021) 12:1337 | https://doi.org/10.1038/s41467-021-21583-9 | www.nature.com/naturecommunications,](https://www.nature.com/articles/s41467-021-21583-9)   
*


Documentation
* [SnapATAC2 Github page](https://github.com/kaizhang/SnapATAC2)
* [SnapATAC2 Tutiorials](https://kzhang.org/SnapATAC2/tutorials/index.html)


Important Notes
* Module Name: SnapATAC2 (see [the modules page](/apps/modules.html) for more information)
* Unusual environment variables set
	+ **SNAPATAC2\_HOME**   installation directory
	+ **SNAPATAC2\_BIN**         executable directory
	+ **SNAPATAC2\_SRC**     test data directory
	+ **SNAPATAC2\_TESTS**     test scripts folder



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive** 
[user@cn4199 ~]$ **module load SnapATAC2**
[+] Loading singularity  3.10.3  on cn4199
[+] Loading SnapATAC2 2.1.2  ...

```

Running the test scripts:

```

[user@cn4199 ~]$ **python $SNAPATAC2\_TESTS/test\_func.py**
[user@cn4199 ~]$ **python $SNAPATAC2\_TESTS/test\_pickle.py**
[user@cn4199 ~]$ **python $SNAPATAC2\_TESTS/test\_similarity.py**
[user@cn4199 ~]$ **python $SNAPATAC2\_TESTS/test\_tools.py**

```

Importing the python module snapatac2:

```

[user@cn4199 ~]$ **python**
Python 3.7.12 | packaged by conda-forge | (default, Oct 26 2021, 06:08:53)
[GCC 9.4.0] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> import snapatac2 as snap
>>> snap.__version__
'2.1.2'

```





