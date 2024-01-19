

document.querySelector('title').textContent = 'FEBio: Finite Elements for Biomechanics';
FEBio: Finite Elements for Biomechanics


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
[Swarm of jobs](#swarm) 
 |



FEBio software suite implement a nonlinear implicit finite element (FE) framework, 
designed specifically for analysis in computational solid biomechanics. 
FEBio offers modeling scenarios, constitutive models, and boundary conditions, 
which are relevant to numerous applications in biomechanics. 
The open-source FEBio software is written in C++, 
with particular attention to scalar and parallel performance on modern computer architectures. 



### References:


* Steve A. Maas, Benjamin J. Ellis, Gerard A. Ateshian and Jeffrey A. Weiss   

*FEBio: Finite Elements for Biomechanics*    

[Journal of Biomechanical Engineering, Vol.134 (2012)](https://asmedigitalcollection.asme.org/biomechanical/article/134/1/011005/455684/FEBio-Finite-Elements-for-Biomechanics)


Documentation
* [FEBio home page](http://mrl.sci.utah.edu/software/febio)


Important Notes
* Module Name: FEBio (see [the modules page](/apps/modules.html) for more information)
* Unusual environment variables set
	+ **FEBIO\_HOME**  installation directory
	+ **FEBIO\_BIN**  executable directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive -c 6 --mem 15g --gres=lscratch:20**
[user@cn3144 ~]$ **module load FEBio**
[+] Loading FEBio  v4.0
[user@cn3144 ~]$  **febio4** 
===========================================================================
         ________    _________   _______       __     _________
        |        |\ |        |\ |       \\    |  |\  /         \\
        |    ____|| |    ____|| |    __  ||   |__|| |    ___    ||
        |   |\___\| |   |\___\| |   |\_| ||    \_\| |   //  \   ||
        |   ||__    |   ||__    |   ||_| ||   |  |\ |  ||    |  ||
        |       |\  |       |\  |         \\  |  || |  ||    |  ||
        |    ___||  |    ___||  |    ___   || |  || |  ||    |  ||
        |   |\__\|  |   |\__\|  |   |\__|  || |  || |  ||    |  ||
        |   ||      |   ||___   |   ||__|  || |  || |   \\__/   ||
        |   ||      |        |\ |          || |  || |           ||
        |___||      |________|| |_________//  |__||  \_________//

      F I N I T E   E L E M E N T S   F O R   B I O M E C H A N I C S

  version 4.0.0
  FEBio is a registered trademark.
  copyright (c) 2006-2022 - All rights reserved

===========================================================================

Default linear solver: skyline

febio>
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```





