

document.querySelector('title').textContent = 'ITK-SNAP-BS: segmentation of structures in 3D and 4D biomedical images.';
**ITK-SNAP-BS: segmentation of structures in 3D and 4D biomedical images.** 


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



ITK-SNAP brings active contour
segmentation to the fingertips of clinical researchers.
This application fulfills a specific and pressing need of biomedical
imaging research by providing a combination of manual and
semiautomatic tools for extracting structures in 3D image data of
different modalities and from different anatomical regions.



### References:


* Paul A. Yushkevich, Joseph Piven, Heather Cody Hazlett, Rachel Gimpel Smith, Sean Ho,
James C. Gee and Guido Gerig   

*User-guided 3D active contour segmentation of anatomical structures:
Significantly improved efficiency and reliability.*   

[NeuroImage 31 (2006) 1116 – 1128](https://www.sciencedirect.com/science/article/pii/S1053811906000632)


Documentation
* [ITK-SNAP-Burgess-Sodt Github page](https://github.com/alexsodt/ITK-SNAP-Burgess-Sodt)
* [ITK-snap-builder Github page](https://github.com/bsandbro/ITK-snap-builder/)
* [ITK-SNAP Documentation](http://www.itksnap.org/pmwiki/pmwiki.php?n=Documentation.SNAP3)


Important Notes
* Module Name: ITK-SNAP-BS (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set
	+ **IS\_HOME**  home folder
	+ **IS\_BIN**   executable folder
	+ **IS\_SRC**   source folder



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
  

 As the first step, please allocate an [interactive session](/docs/userguide.html#int):


The use of ITK-SNAP-BS GUI requires a 
[graphical X11 connection](https://hpc.nih.gov/docs/connect.html).   

Both NX and MobaXterm work well for Windows users,   
 while XQuartz works well for Mac users.



```

[user@biowulf]$ **sinteractive --gres=gpu:p100:1,lscratch:10 --mem=20g -c14**
[user@cn4469 ~]$ **module load DeepLabCut/2.2.2** 
[+] Loading singularity  3.10.5  on cn0825
[+] Loading ITK-SNAP-BS  20131007

```

Invoke the ITK-SNAP-BS GUI:

```

[user@cn4469 user]$ **itksnap**

```

![](ITK-SNAP-BS/GUI.png)
  
  
...





