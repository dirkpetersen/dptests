

document.querySelector('title').textContent = 'straw: rapidly stream data from .hic files';
straw: rapidly stream data from .hic files


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



Straw is library which allows rapid streaming of contact data from .hic files.
  




### References:


* Neva C. Durand, James T. Robinson, Muhammad S. Shamim, Ido Machol, Jill P. Mesirov,
Eric S. Lander, and Erez Lieberman Aiden1   

*Juicebox Provides a Visualization System for Hi-C Contact Maps with Unlimited Zoom.*    

[Cell Systems **3** (2016), 99–101](https://www.sciencedirect.com/science/article/pii/S240547121500054X)
DOI: 10.1016/j.cels.2015.07.012


Documentation
* [straw Github page](https://github.com/aidenlab/straw)


Important Notes
* Module Name: straw (see [the modules page](/apps/modules.html) for more information)
* Unusual environment variables set
	+ **STRAW\_HOME**  installation directory
	+ **STRAW\_BIN**  executable directory
	+ **STRAW\_DATA**  sample data directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive**
[user@cn3144 ~]$ **module load straw**
[+] Loading singularity  3.10.5  on cn4172
[+] Loading straw  1.3.1

```

  

The **straw** application installed on Biowulf comprises the R and Python versions.   

  

The basdic usage of the R version is as follows:

```

[user@cn3144 ~]$ **R**
R version 4.2.3 (2023-03-15) -- "Shortstop Beagle"
...
> library(strawr)
> library(mariner)
>
...

```

  

To use the Python version, type:

```

[user@cn3144 ~]$ **python**
Python 3.10.9 (main, Mar  8 2023, 10:47:38) [GCC 11.2.0] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> import hicstraw
>>> import mustache
>>>
...

```

  

End the interactive session:

```

[user@cn3111 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```





