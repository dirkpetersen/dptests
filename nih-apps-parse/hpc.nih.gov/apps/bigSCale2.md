

document.querySelector('title').textContent = 'bigSCale2: Framework for clustering, phenotyping, pseudotiming and inferring gene regulatory networks from single cell data';
**bigSCale2: Framework for clustering, phenotyping, pseudotiming and inferring gene regulatory networks from single cell data**


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



bigSCale is a complete framework for the analysis and visualization of single cell data. It allows to cluster, phenotype, perform pseudotime analysis, infer gene regulatory networks and reduce large datasets in smaller datasets with higher quality.



### References:


* Giovanni Iacono, Elisabetta Mereu, Amy Guillaumet-Adkins, Roser Corominas, Ivon Cuscó, Gustavo Rodríguez-Esteban, Marta Gut, Luis Alberto Pérez-Jurado, Ivo Gut and Holger Heyn  

*bigSCale: an analytical framework for big-scale single-cell data*   

[Genome Research](https://genome.cshlp.org/content/28/6/878.full) **28**:878–890.


Documentation
* [bigSCale GitHub page](https://github.com/dfajar2/bigSCale)
* [bigSCale2 GitHub page](https://github.com/iaconogi/bigSCale2)


Important Notes
* When you execute the bigSCale2 command, a new shell is opened for you within a Singularity container. Your environment will change and you will have access to a different set of commands and executables. Please remember to exit this new shell when you are finished with your session.
* Module Name: bigSCale2 (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Environment variables set
	+ **BIGSCALE2\_HOME**  bigSCale2 installation directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf ~]$ **sinteractive --mem=4g**
salloc.exe: Pending job allocation 49998864
salloc.exe: job 49998864 queued and waiting for resources
salloc.exe: job 49998864 has been allocated resources
salloc.exe: Granted job allocation 49998864
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn0868 are ready for job
srun: error: x11: no local DISPLAY defined, skipping

[user@cn0868 ~]$ **module load bigSCale2**
[+] Loading singularity  3.5.3  on cn0868
[+] Loading bigSCale2 20191119  ...

[user@cn0868 ~]$ **bigSCale2**

INFO: You are now opening a new shell within a Singularity container. Your
INFO: environment will change and you will have access to a different set of
INFO: commands and executables. Please remember to exit this shell when you are
INFO: finished with your session.

WARNING: Bind mount '/home/user => /home/user' overlaps container CWD /home/user, may not be available
Singularity> **R**

R version 3.6.1 (2019-07-05) -- "Action of the Toes"
Copyright (C) 2019 The R Foundation for Statistical Computing
Platform: x86_64-pc-linux-gnu (64-bit)

R is free software and comes with ABSOLUTELY NO WARRANTY.
You are welcome to redistribute it under certain conditions.
Type 'license()' or 'licence()' for distribution details.

R is a collaborative project with many contributors.
Type 'contributors()' for more information and
'citation()' on how to cite R or R packages in publications.

Type 'demo()' for some demos, 'help()' for on-line help, or
'help.start()' for an HTML browser interface to help.
Type 'q()' to quit R.

> **library(bigSCale)**

[snip...]

> **quit()**
Save workspace image? [y/n/c]: **n**
Singularity>  exit
exit

[user@cn0868 ~]$ **exit**
exit
salloc.exe: Relinquishing job allocation 49998864

[user@biowulf ~]$ 

```





