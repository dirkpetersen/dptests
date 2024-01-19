

document.querySelector('title').textContent = 'ctffind on Biowulf';
ctffind on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
 |



CTFFIND and CTFTILT are two programs for finding CTFs of electron micrographs.



### References:


* [Mindell, JA, Grigorieff N. 2003. **Accurate determination of local defocus and specimen tilt in electron microscopy.** *J Struct Biol. 142:334-47.*](https://www.ncbi.nlm.nih.gov/pubmed/12781660)
* [Rohou, A, Grigorieff N. 2015. **CTFFIND4: Fast and accurate defocus estimation from electron micrographs.** *J Struct Biol. 192:216–221.*](https://www.ncbi.nlm.nih.gov/pubmed/26278980)


Documentation
* [CTFFIND Main Page (Grigorieff Lab)](http://grigoriefflab.janelia.org/ctf)


Important Notes
* Module Name: ctffind (see [the modules page](/apps/modules.html) for more information)
* Multithreaded



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load ctffind**
[user@cn3144 ~]$ **ctffind --help**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. ctffind.sh). For example:



```
#!/bin/bash
module load ctffind

time ctffind3.exe << eof
micrograph.mrc
montage.pow
2.0,200.0,0.07,60000,7.0                        !CS[mm],HT[kV],AmpCnst,XMAG,DStep[um]
128,200.0,8.0,5000.0,30000.0,1000.0,100.0       !Box,ResMin[A],ResMax[A],dFMin[A],dFMax[A],FStep[A],dAst[A]
eof

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] ctffind.sh
```







