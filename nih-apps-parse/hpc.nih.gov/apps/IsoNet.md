

document.querySelector('title').textContent = 'RGI: Resistance Gene Identifier ';
**IsoNet: Isotropic Reconstruction of Electron Tomograms with Deep Learning**


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


  


IsoNet is a deep learning-based software package that iteratively reconstructs the missing-wedge information and increases signal-to-noise ratio, using the knowledge learned from raw tomograms. Without the need for sub-tomogram averaging, IsoNet generates tomograms with significantly reduced resolution anisotropy. Applications of IsoNet to three representative types of cryoET data demonstrate greatly improved structural interpretability: resolving lattice defects in immature HIV particles, establishing architecture of the paraflagellar rod in Eukaryotic flagella, and identifying heptagon-containing clathrin cages inside a neuronal synapse of cultured cells. detect AMR genes from thirteen genomes of Pseudomonas strains.



  

### References:


* Yun-Tao Liu, Heng Zhang, Hui Wang, Chang-Lu Tao, Guo-Qiang Bi & Z. Hong Zhou   

*Isotropic reconstruction for electron tomography with deep learning*  

[Nature Communications volume 13, Article number: 6482 (2022)](https://www.nature.com/articles/s41467-022-33957-8)


Documentation
* [IsoNet on Github](https://github.com/IsoNet-cryoET/IsoNet)


Important Notes
* Module Name: IsoNet (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set
	+ **ISONET\_HOME**  installation directory
	+ **ISONET\_BIN**       bin directory
	+ **ISONET\_SRC**       source directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive** 
[user@cn3107 ~]$ **module load IsoNet**
[+] Loading singularity  3.10.5  on cn4183
[+] Loading CUDA Toolkit  10.2.89  ...
[+] Loading cuDNN/7.6.5/CUDA-10.2 libraries...
[+] Loading IsoNet  0.2.1
[user@cn3107 ~]$ **isonet.py -h** 
INFO: Showing help with the command 'isonet.py -- --help'.

NAME
    isonet.py - ISONET: Train on tomograms and restore missing-wedge

SYNOPSIS
    isonet.py -

DESCRIPTION
    for detail discription, run one of the following commands:

    isonet.py prepare_star -h
    isonet.py prepare_subtomo_star -h
    isonet.py deconv -h
    isonet.py make_mask -h
    isonet.py extract -h
    isonet.py refine -h
    isonet.py predict -h
    isonet.py resize -h
    isonet.py gui -h

[user@cn3107 ~]$ **isonet.py gui** 
![](isonet/isonet_gui.png)

```

End the interactive session:

```

[user@cn3107 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```





