

document.querySelector('title').textContent = 'NovaCTF on Biowulf';
NovaCTF on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
 |



Description



### References:


* Turoňová, B., Schur, F.K.M, Wan, W. and Briggs, J.A.G.
 [**Efficient 3D-CTF correction for cryo-electron tomography using NovaCTF improves subtomogram averaging resolution to 3.4 Å.**](https://doi.org/10.1016/j.jsb.2017.07.007)
*J Struct Biol. 2017 Sep;199(3):1870-195*


Documentation
* [NovaCTF Main Site](https://github.com/turonova/novaCTF)


Important Notes
* Module Name: novactf (see [the modules page](/apps/modules.html) for more information)
 * singlethreaded
 * Environment variables set 
	+ $NOVACTF\_HOME
	+ Example files in $NOVACTF\_HOME/setup\_examples


This app requires the clip command from [IMOD](/apps/IMOD.html).


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

[user@cn3144 ~]$ **module load novactf IMOD**

[user@cn3144 ~]$ **clip flipyz input.ali output.ali**

[user@cn3144 ~]$ **novaCTF -Algorithm defocus -InputProjections input\_stack.st -FULLIMAGE 464,464 -THICKNESS 140 -TILTFILE angles.tlt -SHIFT 0.0,0.0 -CorrectionType phaseflip -DefocusFileFormat ctffind4 -CorrectAstigmatism 1 -DefocusFile defocus\_file.txt -PixelSize 0.135 -DefocusStep 15 -DefocusShiftFile file\_with\_additional\_defocus.txt**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. novactf.sh). For example:



```

#!/bin/bash
set -e
module load novactf IMOD
novaCTF -Algorithm defocus -InputProjections input_stack.st -FULLIMAGE 464,464 -THICKNESS 140 -TILTFILE angles.tlt -SHIFT 0.0,0.0 -CorrectionType phaseflip -DefocusFileFormat ctffind4 -CorrectAstigmatism 1 -DefocusFile defocus_file.txt -PixelSize 0.135 -DefocusStep 15 -DefocusShiftFile file_with_additional_defocus.txt

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] novactf.sh
```





