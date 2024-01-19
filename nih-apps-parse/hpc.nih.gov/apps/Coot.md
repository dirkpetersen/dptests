

document.querySelector('title').textContent = 'Coot on Biowulf';
Coot on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
 |


Coot is for macromolecular model building, model completion and validation, particularly suitable for protein modelling using X-ray data.


Coot displays maps and models and allows model manipulations such as idealization, real space refinement, manual rotation/translation, rigid-body fitting, ligand search, solvation, mutations, rotamers, Ramachandran plots, skeletonization, non-crystallographic symmetry and more.


### References:


* [Emsley P, Lohkamp B, Scott WG, Cowtan K. Features and development of Coot. *Acta Crystallogr D Biol Crystallogr. 2010 Apr* 66(Pt 4):486-501.](https://www.ncbi.nlm.nih.gov/pubmed/?term=20383002)


Documentation
* [Coot documentation](http://www2.mrc-lmb.cam.ac.uk/personal/pemsley/coot/web/docs/coot.html)


Important Notes
* Module Name: Coot (see [the modules page](/apps/modules.html) for more information)
 * Singlethreaded


This application requires an [X-Windows connection](/docs/connect.html).


Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load Coot**
[user@cn3144 ~]$ **coot**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```



