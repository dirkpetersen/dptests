

document.querySelector('title').textContent = 'PyTom: a toolbox developed for interpreting cryo electron tomography data';
**PyTom: a toolbox developed for interpreting cryo electron tomography data**


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



PyTom is a software package for the analysis of volumetric data 
obtained by cryo electron tomography (cryo-ET). 
It covers a complete pipeline of processing steps 
for tomogram reconstruction, localization of macromolecular complexes in tomograms, 
fine alignment of subtomograms extracted at these locations, and their classification.



### References:


* Thomas Hrabe, Yuxiang Chen, Stefan Pfeffer, Luis Kuhn Cuellar, Ann-Victoria Mangold   
& Friedrich Förster   

 *PyTom: A python-based toolbox for localization of macromolecules in cryo-electron
tomograms and subtomogram analysis*   

[Journal of Structural Biology 178 (2012) 177–188](https://reader.elsevier.com/reader/sd/pii/S1047847711003492?token=C0FF8D1654C70F8ABC07CFCE18C3E62AC035C09981B31F6197FA7094484F05BC5BB5D70DC5A4CDB95CEB2F7F4B1F9AEA&originRegion=us-east-1&originCreation=20221018175925)


Documentation
* [PyTom Github page](https://github.com/FridoF/PyTom)


Important Notes
* Module Name: PyTom (see [the modules page](/apps/modules.html) for more information)
* Unusual environment variables set
	+ **PYTOM\_HOME**   installation directory
	+ **PYTOM\_BIN**         executable directory
	+ **PYTOM\_DATA**     sample data directory
	+ **PYTOM\_TESTS**     test scripts folder



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive** 
[user@cn4199 ~]$ **module load PyTom**
[+] Loading PyTom 1.0  ...

```

Running tests:

```

[user@cn4199 ~]$ **wget https://github.com/FridoF/PyTom/archive/refs/tags/v1.0.tar.gz**
[user@cn4199 ~]$ **tar -zxf v1.0.tar.gz && rm -f v1.0.tar.gz**
[user@cn4199 ~]$ **cd PyTom-1.0/tests**
[user@cn4199 ~]$ **pytom -m unittest discover**
...
This license affects the software package PyTom and all the herein distributed source / data files.
Authors:
Marten Chaillet
Gijs van der Schot
Ilja Gubins
Mihajlo Vanevic
Thomas Hrabe
Yuxiang Chen
Friedrich Foerster

Copyright (c) 2021
Utrecht University
http://www.pytom.org

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

The complete license can be obtained from

http://www.gnu.org/licenses/gpl-2.0.html.

RMSD =  0.72
....................We adjusted msdz to 5.000 nm to make it an integer multiple of pixel size.
Number of slices for multislice:  12
Determining DQE for K2SUMMIT
chi-square value for fitting sinc function is 0.0005657000438432647
Determining MTF for K2SUMMIT
chi-square value for fitting sinc function is 1.0941206081582777e-05
Projecting the model with 1 processes
[Parallel(n_jobs=1)]: Using backend SequentialBackend with 1 concurrent workers.
Transforming sample for tilt/frame  0
Reduced rotation height for relevant specimens: 300
Simulating projection with multislice method for frame/tilt  0
Number of electrons per pixel (before oversampling and sample absorption): 100.00
[Parallel(n_jobs=1)]: Done   1 out of   1 | elapsed:    8.2s remaining:    0.0s
[Parallel(n_jobs=1)]: Done   1 out of   1 | elapsed:    8.2s finished
...
.variance for spline interpolation with numba: 3.850465873256326e-06
.before parallel:  24
after setting for parallel:  4
execution numba 4 threads:  0.1669407606124878
variance for spline interpolation with numba parallel: 3.816609478235478e-06
.variance for pt linear interpolation: 3.833221927715919e-06
variance for pt cubic interpolation: 3.817091292496215e-06
variance for pt spline interpolation: 3.816610180137391e-06
.variance for vt linear interpolation on cpu: 3.8332236727001145e-06
variance for vt bspline interpolation on cpu: 3.826985903288005e-06
variance for vt bspline_simple interpolation on cpu: 3.826985903288005e-06
variance for vt filt_bspline interpolation on cpu: 3.910493887815392e-06
variance for vt filt_bspline_simple interpolation on cpu: 3.910493887815392e-06
.variance for vt linear interpolation on gpu: 3.846456e-06
variance for vt bspline interpolation on gpu: 3.8316716e-06
variance for vt bspline_simple interpolation on gpu: 3.8315306e-06
variance for vt filt_bspline interpolation on gpu: 3.879772e-06
variance for vt filt_bspline_simple interpolation on gpu: 3.87833e-06
...
----------------------------------------------------------------------
Ran 55 tests in 72.486s

OK

[user@cn4199 ~]$ **exit**
salloc.exe: Relinquishing job allocation 59748321
[user@biowulf ~]$

```





