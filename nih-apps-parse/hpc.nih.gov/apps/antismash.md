

document.querySelector('title').textContent = 'antiSMASH on Biowulf';
antiSMASH on Biowulf


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



 antiSMASH allows the rapid genome-wide identification, annotation and analysis of secondary metabolite biosynthesis gene clusters in bacterial and fungal genomes.
 It integrates and cross-links with a large number of in silico secondary metabolite analysis tools that have been published earlier.



### References:


* antiSMASH: Rapid identification, annotation and analysis of secondary metabolite biosynthesis gene clusters.
 Marnix H. Medema, Kai Blin, Peter Cimermancic, Victor de Jager, Piotr Zakrzewski, Michael A. Fischbach, Tilmann Weber, Rainer Breitling & Eriko Takano
 Nucleic Acids Research (2011) doi: [10.1093/nar/gkr466](https://doi.org/10.1093/nar/gkr466)
* antiSMASH 6.0: improving cluster detection and comparison capabilities.
 Kai Blin, Simon Shaw, Alexander M Kloosterman, Zach Charlop-Powers, Gilles P van Weezel, Marnix H Medema, & Tilmann Weber
 Nucleic Acids Research (2021) doi: [10.1093/nar/gkab335](https://doi.org/10.1093/nar/gkab335).


Documentation
* [antiSMASH Documentation](https://docs.antismash.secondarymetabolites.org)
* [antiSMASH on GitHub](https://github.com/antismash/antismash)


Important Notes
* Module Name: antismash (see [the modules page](/apps/modules.html) for more information)
 * Multithreaded via the -c/--cpus argument.
 As of version 6.1, the default is **128 cpus**, so make sure to set this parameter to what you have allocated.
 * This application produces HTML reports. You can use [hpcdrive](/docs/hpcdrive.html) to view these reports on your local workstation.
* Environment variables set
