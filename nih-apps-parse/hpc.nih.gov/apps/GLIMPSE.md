

document.querySelector('title').textContent = 'GLIMPSE: Genotype Likelihoods IMputation and PhaSing mEthod ';
**GLIMPSE: Genotype Likelihoods IMputation and PhaSing mEthod** 


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



GLIMPSE is a phasing and imputation method for large-scale low-coverage sequencing studies.
It perform accurate imputed genotype calls and outperforms SNP arrays.



### References:


* S. Rubinacci, D.M. Ribeiro, R. Hofmeister, O. Delaneau   

*Efficient phasing and imputation of low-coverage sequencing data using large reference panels*   

[bioRxiv preprint doi: https://doi.org/10.1101/2020.04.14.040329;](https://www.biorxiv.org/content/10.1101/2020.04.14.040329v1)


Documentation
* [GLIMPSE Home page](https://odelaneau.github.io/GLIMPSE/index.html)


Important Notes
* Module Name: GLIMPSE (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set
	+ **GLIMPSE\_HOME**  installation directory
	+ **GLIMPSE\_BIN**       executable directory
	+ **GLIMPSE\_DATA**  sample data directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=16g --cpus-per-task=16**
[user@cn3335 ~]$**module load GLIMPSE**
[+] Loading GLIMPSE  1.1.1

```


```

[user@biowulf]$ **GLIMPSE\_chunk** 

[GLIMPSE] Split chromosomes into chunks
  * Author        : Simone RUBINACCI & Olivier DELANEAU, University of Lausanne
  * Contact       : simone.rubinacci@unil.ch & olivier.delaneau@unil.ch
  * Version       : 1.1.1
  * Run date      : 21/11/2022 - 17:25:19
...
[user@biowulf]$ **GLIMPSE\_ligate** 
...
[user@biowulf]$ **GLIMPSE\_sample**
...
[user@biowulf]$ **GLIMPSE\_concordance**
...
[user@biowulf]$ **GLIMPSE\_phase**
...
[user@biowulf]$ **GLIMPSE\_snparray**
...

```

[user@cn3335 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$





