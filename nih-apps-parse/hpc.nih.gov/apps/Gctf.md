

document.querySelector('title').textContent = 'Gctf on Biowulf';
Gctf on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
 |



Gctf provides accurate estimation of the contrast transfer function (CTF) for near-atomic resolution cryo electron microscopy (cryoEM) reconstruction using GPUs. The main target of Gctf is to maximize the cross-correlation of a simulated CTF with the logarithmic amplitude spectra (LAS) of observed micrographs after background subtraction.



### References:


* Zhang K.
[**Gctf: Real-time CTF determination and correction.**](https://www.ncbi.nlm.nih.gov/pubmed/26592709)
*J Struct Biol. 2016 Jan;193(1):1-12.*


Documentation
* [Gctf Main Site](https://www.mrc-lmb.cam.ac.uk/kzhang/)
* Type **Gctf --help**


Important Notes
* Module Name: Gctf (see [the modules page](/apps/modules.html) for more information)
* GPU-accelerated
* environment variables set 
	+ GCTF\_HOME
	+ RELION\_GCTF\_EXECUTABLE


Gctf can utilize GPUs. This requires the user load the proper CUDA library for the executables. For example, the executable **Gctf-v1.06\_sm\_30\_cu7.5\_x86\_64** requires **module load CUDA/7.5** prior to use.


Alternatively, a wrapper script which automatically loads the correct CUDA library can be used: **Gctf**


Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int), along with at least one GPU, and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --constraint=gpup100 --gres=gpu:p100:1**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load Gctf CUDA/7.5**
[user@cn3144 ~]$ Gctf-v1.06_sm_30_cu7.5_x86_64 --apix 1.07  --kV 300 --Cs 2.7 --ac 0.1  Micrographs/Falcon*.mrc

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. Gctf.sh). For example:



```

#!/bin/bash
module load Gctf
Gctf --apix 1.07 --kV 300 --Cs 2.7 --ac 0.1 Micrographs/Falcon*.mrc

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] --gres=gpu:p100:1 --partition=gpu Gctf.sh
```





