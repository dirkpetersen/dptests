

document.querySelector('title').textContent = 'CRYOLO on Biowulf';
CRYOLO on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
 |



CRYOLO is a fast and accurate particle picking procedure based on convolutional neural networks.



### References:


* [SPHIRE-crYOLO is a fast and accurate fully automated particle picker for cryo-EM.](https://doi.org/10.1038/s42003-019-0437-z) Communications Biology 2, Wagner, T. et al. (2019).


Documentation
* [cryolo Main Site](http://sphire.mpg.de/wiki/doku.php?id=pipeline:window:cryolo)


Important Notes
* Module Name: cryolo (see [the modules page](/apps/modules.html) for more information)



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

[user@cn3144 ~]$ **module load cryolo**
[+] Loading cryolo  1.5.3 

[user@cn3144 ~]$ **cryolo\_gui.py &**

```

The above should display the Graphical User Interface of cryolo







