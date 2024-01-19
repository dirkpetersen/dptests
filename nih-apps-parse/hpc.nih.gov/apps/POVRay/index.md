

document.querySelector('title').textContent = 'POVRay on Biowulf';
POVRay on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
 |



POV-Ray is a high-quality tool for creating three-dimensional graphics.
Raytraced images are publication-quality and 'photo-realistic', but are
computationally expensive so that large images can take many hours to create.
POV-Ray images can also require more memory than many desktop machines can
handle. To address these concerns, a parallelized version of POV-Ray has been
installed on Biowulf.



POV-Ray output is limited to only .png, .tga, or .ppm image formats. There
are number of programs which can convert images from one format to another
available on Helix Systems (e.g., convert, [gimp](http://www.gimp.org), [imagemagick](http://www.imagemagick.org/script/index.php), [xnview](http://www.xnview.com)).


Documentation
* [POV-Ray Main Site](http://www.povray.org/)
* Type **man povray** after loading the module


Important Notes
* Module Name: POVRay (see [the modules page](/apps/modules.html) for more information)
* Multithreaded/singlethreaded/MPI...
* Environment variables set 
	+ PATH
	+ MANPATH
	+ POVRAY\_EXAMPLES* Example files in $POVRAY\_EXAMPLES


By default, POV-Ray attempts to display the image as it is being rendered to the screen. This feature requires an [X-Windows connection](/docs/connect.html). To disable the display, you must include **-D** in the POV-Ray command line.


Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive --cpus-per-task=16 --threads-per-core=1**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **cp $POVRAY\_EXAMPLES/1asy.pov .**
[user@cn3144 ~]$ **povray +H542 +W456 -I1asy.pov -O1asy.pov.tga +P +X +A +FT +C -D +wt${SLURM\_CPUS\_PER\_TASK}**
...
Render Time:
  Photon Time:      No photons
  Radiosity Time:   No radiosity
  Trace Time:       0 hours  0 minutes  0 seconds (0.349 seconds)
              using 30 thread(s) with 8.276 CPU-seconds total
POV-Ray finished

```

Once the image render is completed, convert the file into a JPEG and then display it via X11:



```

[user@cn3144 ~]$ **convert 1asy.pov.tga 1asy.jpg**
[user@cn3144 ~]$ **display 1asy.jpg**
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```

![1asy image](1asy.jpg)



Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. POVRay.sh). For example:



```

#!/bin/bash
module load POVRay
povray +H1125 +W858 -Iag9.pov -O1ag9.pov.tga +P +X +A +FT +C +wt${SLURM_CPUS_PER_TASK}
povray +H2170 +W1826 -I1asy.pov -O1asy.pov.tga +P +A +FT +wt${SLURM_CPUS_PER_TASK}

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command. **--threads-per-core=1** disables hyperthreading on the node.



```
sbatch [--cpus-per-task=#] [--mem=#] --threads-per-core=1 POVRay.sh
```





