

document.querySelector('title').textContent = 'PyMOL on Biowulf';
PyMOL on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[PoseFilter](#posefilter)
[Interactive job](#int) 
 |



Description
PyMOL is a powerful and comprehensive molecular visualization product for rendering and animating 3D molecular structures. Pymol is a user-sponsored molecular visualization system on an open-source foundation.

Documentation
* PyMOL Main Site: [PyMOL on Sourceforge](https://sourceforge.net/projects/pymol)


Important Notes
This application requires a [graphical connection using NX](/docs/connect.html#nx)


* Module Name: pymol (see [the modules page](/apps/modules.html) for more information)
 * Must have NX software installed: [NX](/docs/nx.html)* Must use an interactive node
 * To open pdb files or fetch from pdb, turn shaders off:
	+ Start pymol
	 + Type: set use\_shaders, 0
	 + Type: fetch xxxx* The shaders setting can be made permanent by editing ~/.pymolrc:  


```
$ echo 'set use_shaders, 0' >> ~/.pymolrc
```


PyMOL will **crash** unless shaders are turned off.



PoseFilter
There is a special install of PyMOL that includes the [PoseFilter plugin](https://github.com/skalyaanamoorthy/PoseFilter). It can be used by loading the posefilter module.


Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --cpus-per-tasks=8 --mem=#g**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load pymol**
[user@cn3144 ~]$ **pymol&**
[user@cn3144 ~]$ libGL error: unable to load driver: swrast_dri.so
libGL error: failed to load driver: swrast
 PyMOL(TM) Molecular Graphics System, Version 1.8.6.0.
 Copyright (c) Schrodinger, LLC.
 All Rights Reserved.
 
    Created by Warren L. DeLano, Ph.D. 
 
    PyMOL is user-supported open-source software.  Although some versions
    are freely available, PyMOL is not in the public domain.
 
    If PyMOL is helpful in your work or study, then please volunteer 
    support for our ongoing efforts to create open and affordable scientific
    software by purchasing a PyMOL Maintenance and/or Support subscription.

    More information can be found at "http://www.pymol.org".
 
    Enter "help" for a list of commands.
    Enter "help " for information on a specific command.

 Hit ESC anytime to toggle between text and graphics.

 Detected OpenGL version prior to 2.0. Shaders and volumes unavailable.
 OpenGL graphics engine:
 GL\_VENDOR: Mesa Project
 GL\_RENDERER: Software Rasterizer
 GL\_VERSION: 1.4 (2.1 Mesa 10.0.1)
 Detected 8 CPU cores. Enabled multithreaded rendering.

```

You should see the PyMOL GUI windows appear on your screen.
![Pymol display](PyMOL.png)












