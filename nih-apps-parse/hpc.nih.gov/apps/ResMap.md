

document.querySelector('title').textContent = 'ResMap on Biowulf';
ResMap on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int)
[Batch job](#sbatch)
[GPU Acceleration](#gpus)
 |


ResMap (Resolution Map) is a Python (NumPy/SciPy) application with a Tkinter GUI
and a command-line interface. It is a software package for computing the local resolution
of 3D density maps studied in structural biology, primarily electron cryo-microscopy (cryo-
EM).


### References:


* A. Kucukelbir, F.J. Sigworth, H.D. Tagare.
 [**Quantifying the Local Resolution of Cryo-EM Density Maps.**](http://www.ncbi.nlm.nih.gov/pubmed/24213166)
*Nature Methods, Volume 11, Issue 1, Pages 63-65, 2014.*


Documentation
* ResMap Main Site: <http://resmap.sourceforge.net/>
* ResMap PDF: <ResMap-manual_v-1.95.pdf>


Important Notes
* Module Name: ResMap (see [the modules page](/apps/modules.html) for more information)


This application requires an [X-Windows connection](/docs/connect.html). Users are encouraged to use [NX](https://hpc.nih.gov/docs/nx.html) as their X11 servers.


Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --gres=gpu:p100:1**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load ResMap**

[user@cn3144 ~]$ **ResMap**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```

After typing 'ResMap' at the prompt, input is menu-driven.


![ResMap window](ResMap.PNG)
The correct path for the CUDA library file is **/opt/resmap/ResMap\_krnl-cuda-V9.0.102-sm35\_gpu.so**. You will need to update the GUI.



Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. ResMap.sh). ResMap can be run without the need for an X11 server by starting a dummy X11 server with Xvfb and setting the $DISPLAY variable to a bogus number.



```
#!/bin/bash
module load ResMap
Xvfb -shmem -screen 0 1280x1024x24 &
export DISPLAY=":0"
cp $RESMAP_EXAMPLES/*.map .
ResMap --doBenchMarking --noguiSplit \
  emd_8731_half_map_1.map emd_8731_half_map_2.map
```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --mem=10 ResMap.sh
```

GPU Acceleration
ResMap can be accelerated at least 10-fold by using GPUs rather than CPUs. This requires three options, **--use\_gpu**, **--set\_gpu** and **--lib\_krnl\_gpu**. Here is a batch script for doing so:



```
#!/bin/bash
module load ResMap
Xvfb -shmem -screen 0 1280x1024x24 &
export DISPLAY=":0"
cp $RESMAP_EXAMPLES/*.map .
ResMap --doBenchMarking --noguiSplit \
  emd_8731_half_map_1.map emd_8731_half_map_2.map \
  --use_gpu=yes --set_gpu=0 \
  --lib_krnl_gpu=/opt/resmap/ResMap_krnl-cuda-V9.0.102-sm35_gpu.so
```

Submit this job using this Slurm command, allocating a single GPU for the job:



```
sbatch --mem=10 --gres=gpu:p100:1 ResMap.sh
```





