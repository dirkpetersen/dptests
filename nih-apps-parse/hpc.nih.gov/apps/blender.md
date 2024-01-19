

document.querySelector('title').textContent = 'Blender on Biowulf';
Blender on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Visual partition job](#svis) 
[Batch job](#sbatch) 
[Swarm of jobs](#swarm) 
 |



Blender is the free and open source 3D creation suite. It supports the entirety of the 3D pipeline—modeling, rigging, animation, simulation, rendering, compositing and motion tracking, even video editing and game creation.



Blender on Biowulf is meant for command-line rendering. The .blend file should be created outside of Biowulf, as none of the nodes have support for OpenGL, required for running the GUI.


Documentation
* [Blender Main Site](https://www.blender.org/)


Important Notes
* Module Name: Blender (see [the modules page](/apps/modules.html) for more information)
 * Multithreaded
 * Example files in BLENDER\_EXAMPLES


**NOTE:** Blender running in commandline/batch mode with **--background** may still attempt to connect graphics to an X11 server. To avoid errors, a "fake" X11 server can be started to intercept these connections:



```
Xvfb -shmem -screen 0 1280x1024x24 & ; export DISPLAY=":0" ; blender --background ...
```

Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive --cpus-per-task=8**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load blender**
[user@cn3144 ~]$ **cp $BLENDER\_EXAMPLES/fishy\_cat.blend .**
[user@cn3144 ~]$ **blender -t ${SLURM\_CPUS\_ON\_NODE} -noaudio --background fishy\_cat.blend --render-output run/output --render-frame 1**
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```

Using the visual partition
[Visual partition jobs](/docs/svis.html) should be used when a graphical application requires GPU acceleration and proper OpenGL handling.
Following the directions to start a [visual partition job](/docs/svis.html), open a console window, load a gpu-enabled version of blender, and launch the blender GUI:



```

[user@cn0655 ~]$ **module load blender/2.82\_gpu**
[user@cn0655 ~]$ **vglrun blender**

```

Opening "Edit --> Preferences" should show that blender is using the GPU available on the visual partition node allocated.



Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. blender.sh). For example:



```

#!/bin/bash
module load blender
blender -t ${SLURM_CPUS_ON_NODE} -noaudio --background file.blend --render-output run/output

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] blender.sh
```





