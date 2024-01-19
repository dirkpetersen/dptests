

document.querySelector('title').textContent = 'pairtools: a command-line framework to process sequencing data from a Hi-C experimen ';
**pairtools: a command-line framework to process sequencing data from a Hi-C experimen** 


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



Pairtools is a simple and fast command-line framework to process sequencing data 
from a Hi-C experiment. Pairtools perform various operations on Hi-C pairs and occupy the middle position in a typical Hi-C data processing pipeline. Pairtools aim to be an all-in-one tool for processing Hi-C pairs.



### References:


* Fan Song, Jie Xu, Jesse Dixon and Feng Yue   

[*Analysis of Hi-C Data for Discovery of Structural Variations in Cancer.*](https://link.springer.com/protocol/10.1007/978-1-0716-1390-0_7) Hi-C Data Analysis. Methods and Protocols (2022). pp 143-161


Documentation
* [pairtools Github page](https://github.com/open2c/pairtools)
* [pairtools Documentatiom](https://pairtools.readthedocs.io/en/latest/index.html)


Important Notes
* Module Name: Funannotate (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set
	+ **PT\_HOME**  installation directory
	+ **PT\_BIN**    executables directory
	+ **PT\_DATA**    sample data directory
	+ **PT\_TESTS**    test scripts directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
  


```

[user@biowulf]$ **sinteractive --gres=lscratch:10 --mem=30g -c4**
[user@cn0861 ~]$ **module load pairtools** 
[+] Loading singularity  3.8.5-1  on cn4207
[+] Loading pairtools  0.3.0

```


```

[user@biowulf]$ **python $PT\_TESTS/test\_dedup.py** 

[user@biowulf]$ **python $PT\_TESTS/test\_filterbycov.py** 

[user@biowulf]$ **python $PT\_TESTS/test\_flip.py** 

[user@biowulf]$ **python $PT\_TESTS/test\_headerops.py** 

[user@biowulf]$ **python $PT\_TESTS/test\_markasdup.py** 

[user@biowulf]$ **python $PT\_TESTS/test\_merge.py** 

[user@biowulf]$ **python $PT\_TESTS/test\_parse.py** 

[user@biowulf]$ **python $PT\_TESTS/test\_select.py** 

[user@biowulf]$ **python $PT\_TESTS/test\_sort.py** 

[user@biowulf]$ **python $PT\_TESTS/test\_split.py** 

[user@biowulf]$ **python $PT\_TESTS/test\_stats.py** 

```

etc.


```

[user@cn0861 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226

```





