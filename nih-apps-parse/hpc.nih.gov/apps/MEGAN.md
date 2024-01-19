

document.querySelector('title').textContent = 'MEGAN on Biowulf';
MEGAN on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int)
 |



MEGAN is a MEtaGenome ANalyzer. It can read Blast output. Normally this program would be run on your own desktop, but MEGAN is a memory hog, and some desktop systems may not have enough memory to handle a large dataset.


MEGAN was developed by [D. H. Huson et al at U. Tübingen](https://uni-tuebingen.de/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/informatik/lehrstuehle/algorithms-in-bioinformatics/software/megan6/).




> 
>  MEGAN6 is a comprehensive toolbox for interactively analyzing microbiome data. All the interactive tools you need in one application.
>  * Taxonomic analysis using the NCBI taxonomy or a customized taxonomy such as SILVA
> * Functional analysis using InterPro2GO, SEED, eggNOG or KEGG
> * Bar charts, word clouds, Voronoi tree maps and many other charts
> * PCoA, clustering and networks
> * Supports metadata
> * MEGAN parses many different types of input
> 
> 
> 


### References:


* Daniel H. Huson, Sina Beier, Isabell Flade, Anna Gorska, Mohamed El-Hadidi, Suparna Mitra, Hans-Joachim Ruscheweyh and Rewati Tappu.
 [**MEGAN Community Edition - Interactive Exploration and Analysis of Large-Scale Microbiome Sequencing Data.**](https://doi.org/10.1371/journal.pcbi.1004957)
*PLOS Computational Biology 12(6): e1004957.*



Documentation
* [MEGAN Manual (PDF)](https://software-ab.informatik.uni-tuebingen.de/download/megan6/manual.pdf) at the U. Tübingen site
* [MEGAN Github](https://github.com/husonlab/megan-ce)


Important Notes
This application requires a [graphical connection using NX](/docs/nx.html)


* Module Name: MEGAN (see [the modules page](/apps/modules.html) for more information)
* Based on author suggestions, the maximum memory that MEGAN can use is 32GB. Keep this in mind when allocating memory for you sinteractive session.




Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an interactive node and run MEGAN there. Sample session requesting an allocation of 10 GB of RAM:

```

[user@biowulf ~]$ sinteractive --mem=10g
salloc.exe: Granted job allocation 143361

[user@cn0124 ~]$ module load MEGAN

[user@cn0124 ~]$ MEGAN

![](/images/megan0.png)

![](/images/megan1.png)




[user@cn0124 ~]$ exit
exit
srun: error: cn0124: task 0: Exited with exit code 130
salloc.exe: Relinquishing job allocation 143361
[user@biowulf ~]$



```














