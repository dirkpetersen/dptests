

document.querySelector('title').textContent = 'regenie: whole genome regression modelling of large genome-wide association studies. ';
**regenie: whole genome regression modelling   
 of large genome-wide association studies.** 


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



regenie is a C++ program for whole genome regression modelling   
 
of large genome-wide association studies. It is developed and supported 
by a team of scientists at the Regeneron Genetics Center.
regenie employs the BGEN library.



### References:


* Joelle Mbatchou, Leland Barnard, Joshua Backman, Anthony Marcketta, Jack A. Kosmicki, Andrey Ziyatdinov, Christian Benner, Colm O'Dushlaine, Mathew Barber, Boris Boutkov, Lukas Habegger, Manuel Ferreira, Aris Baras, Jeffrey Reid, Goncalo Abecasis, Evan Maxwell, Jonathan Marchini.   

*Computationally efficient whole genome regression for quantitative and binary traits*  

 [bioRxiv](https://www.biorxiv.org/content/10.1101/2020.06.19.162354v2.abstract) 
(2020),doi: https://doi.org/10.1101/2020.06.19.162354.   
* Band, G. and Marchini, J.   

*BGEN: a binary file format for imputed genotype and haplotype data*  

[bioRxiv](https://www.biorxiv.org/content/10.1101/308296v2) (2018); 
doi: https://doi.org/10.1101/308296


Documentation
* [regenie Github page](https://github.com/rgcgithub/regenie)
* [regenie\_Tutorial](https://rgcgithub.github.io/regenie/options/)
* [BGEN home page](https://enkre.net/cgi-bin/code/bgen/dir?ci=trunk)


Important Notes
* Module Name: regenie (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set
	+ **REGENIE\_HOME**  installation directory
	+ **REGENIE\_BIN**       executable directory
	+ **REGENIE\_SRC**       source code directory
	+ **REGENIE\_DATA**  sample data and checkpoints directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive** 
[user@cn3101 ~]$**module load regenie/3.0.3** 
[+] Loading singularity  3.10.0  on cn3063
[+] Loading regenie  3.0.3

```

The available executables are:

```

[user@cn3101]$ **ls $REGENIE\_BIN** 
bgenix  cat-bgen  edit-bgen  regenie  zstd

```

In particular, the command line options of the executable regenie are as follows: 

```

[user@cn3101]$ **regenie --help**
              |============================|
              |        REGENIE v3.0.3      |
              |============================|

Copyright (c) 2020-2022 Joelle Mbatchou, Andrey Ziyatdinov and Jonathan Marchini.
Distributed under the MIT License.


Usage:
  /regenie/regenie [OPTION...]

  -h, --help      print list of available options
      --helpFull  print list of all available options

 Main options:
      --step INT                specify if fitting null model (=1) or
                                association testing (=2)
      --bed PREFIX              prefix to PLINK .bed/.bim/.fam files
      --pgen PREFIX             prefix to PLINK2 .pgen/.pvar/.psam files
      --bgen FILE               BGEN file
      --sample FILE             sample file corresponding to BGEN file
      --ref-first               use the first allele as the reference for
...

```

To perform training of the predictor network using this executable, copy sample data to the current folder:

```

[user@cn3101]$ **cp $REGENIE\_DATA/\* .**

```

A sample command to run regenie:

```

[user@cn3101]$ **regenie --bgen example.bgen --out my\_output --step 1 --bsize 200 --phenoFile phenotype\_bin.txt**
Start time: Tue Aug 16 13:24:00 2022

              |============================|
              |        REGENIE v3.0.3      |
              |============================|

Copyright (c) 2020-2022 Joelle Mbatchou, Andrey Ziyatdinov and Jonathan Marchini.
Distributed under the MIT License.

Log of output saved in file : my_output.log

Options in effect:
  --bgen example.bgen \
  --out my_output \
  --step 1 \
  --bsize 200 \
  --phenoFile phenotype_bin.txt

Fitting null model
 * bgen             : [example.bgen]
   -summary : bgen file (v1.2 layout, zlib compressed) with 500 named samples and 1000 variants with 8-bit encoding.
   -index bgi file [example.bgen.bgi]
 * phenotypes       : [phenotype_bin.txt] n_pheno = 2
   -keeping and mean-imputing missing observations (done for each trait)
   -number of phenotyped individuals  = 500
 * number of individuals used in analysis = 500
   -residualizing and scaling phenotypes...done (0ms)
 * # threads        : [55]
 * block size       : [200]
 * # blocks         : [5] for 1000 variants
 * # CV folds       : [5]
 * ridge data_l0    : [5 : 0.01 0.25 0.5 0.75 0.99 ]
 * ridge data_l1    : [5 : 0.01 0.25 0.5 0.75 0.99 ]
 * approximate memory usage : 2MB
 * setting memory...done

Chromosome 1
 block [1] : 200 snps  (4ms)
   -residualizing and scaling genotypes...done (3ms)
   -calc working matrices...done (420ms)
   -calc level 0 ridge...done (79ms)
 block [2] : 200 snps  (2ms)
   -residualizing and scaling genotypes...done (1ms)
   -calc working matrices...done (439ms)
   -calc level 0 ridge...done (79ms)
 block [3] : 200 snps  (2ms)
   -residualizing and scaling genotypes...done (1ms)
   -calc working matrices...done (483ms)
   -calc level 0 ridge...done (81ms)
 block [4] : 200 snps  (3ms)
   -residualizing and scaling genotypes...done (1ms)
   -calc working matrices...done (366ms)
   -calc level 0 ridge...done (78ms)
 block [5] : 200 snps  (2ms)
   -residualizing and scaling genotypes...done (1ms)
   -calc working matrices...done (485ms)
   -calc level 0 ridge...done (78ms)

 Level 1 ridge...
   -on phenotype 1 (Y1)...done (0ms)
   -on phenotype 2 (Y2)...done (0ms)

Output
------
phenotype 1 (Y1) :
  0.01  : Rsq = 0.00292408, MSE = 0.995083<- min value
  0.25  : Rsq = 0.00619743, MSE = 0.998022
  0.5   : Rsq = 0.00679147, MSE = 1.00153
  0.75  : Rsq = 0.00753375, MSE = 1.00367
  0.99  : Rsq = 0.00733694, MSE = 1.01373
  * making predictions...writing LOCO predictions...done (9ms)

phenotype 2 (Y2) :
  0.01  : Rsq = 0.012437, MSE = 0.98745<- min value
  0.25  : Rsq = 0.00739346, MSE = 0.997094
  0.5   : Rsq = 0.00612812, MSE = 1.00169
  0.75  : Rsq = 0.00621549, MSE = 1.00343
  0.99  : Rsq = 0.0082828, MSE = 1.00621
  * making predictions...writing LOCO predictions...done (9ms)

List of blup files written to: [my_output_pred.list]

Elapsed time : 2.66076s
End time: Tue Aug 16 13:24:02 2022

```

End the interactive session:

```

[user@cn3101 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```





