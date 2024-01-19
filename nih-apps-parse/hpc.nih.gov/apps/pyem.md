

document.querySelector('title').textContent = 'UCSF pyem on Biowulf';
UCSF pyem on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
 |



UCSF pyem is a collection of Python modules and command-line utilities for electron microscopy of biological samples.



### References:


* Asarnow, D., Palovcak, E., Cheng, Y.
 [**UCSF pyem v0.5.**](https://doi.org/10.5281/zenodo.3576630)
*(Zenodo:2019)*


Documentation
* [UCSF pyem Main Site](https://github.com/asarnow/pyem)


Important Notes
* Module Name: pyem (see [the modules page](/apps/modules.html) for more information)
* singlethreaded



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load pyem**

[user@cn3144 ~]$ **star.py**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```








