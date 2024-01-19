

document.querySelector('title').textContent = 'Protofilament Refinement on Biowulf';
Protofilament Refinement on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
 |



Protofilament Refinement is a software package to further refine microtubule structures, by aligning individual protofilaments rather than full microtubule segements.
Currently, this program assumes that you have already refined your microtubule structures to a moderately high resolution. A full microtubule refinement suite may be further developed in the future.



Documentation
* [Protofilament Refinement Main Site](https://gitlab.com/gedebs371/protofilament-refinement)
* [Tutorial](https://gitlab.com/gedebs371/protofilament-refinement/-/blob/master/Tutorial.pdf)


Important Notes
This application requires a [graphical connection using NX](/docs/nx.html)


* Module Name: pf\_refinement (see [the modules page](/apps/modules.html) for more information)
 * Multithreaded



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive --cpus-per-task=16**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job
[user@cn3144 ~]$ **module load pf\_refinement**
[user@cn3144 ~]$ **cd /data/*user***
[user@cn3144 *user*]$ **mkdir MyNewProject && cd $\_**
[user@cn3144 MyNewProject]$ **pf\_init\_project**
```

This will create a graphical window:


![pf_init_project](pf_refinement_1.png)
After initializing the project in your /data directory and running through the additional steps, the interactive session is ended like so:



```
[user@cn3144 MyNewProject]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```





