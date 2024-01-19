

document.querySelector('title').textContent = 'Octave on Biowulf';
Octave on Biowulf
[![Octave logo](/images/octave.jpg)](http://www.gnu.org/software/octave/)


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job with GUI](#int)
[Interactive job without GUI](#int)
[Batch job](#sbatch) 
[Swarm of jobs](#swarm) 
[Octave Packages](#pkg)
 |



GNU Octave is a
high-level language, primarily intended for numerical computations. It provides
a convenient command line interface for solving linear and nonlinear problems
numerically, and for performing other numerical experiments using a language
that is mostly compatible with Matlab. It may also be used as a batch-oriented
language.


Many, but not all, Matlab .m scripts will run in Octave. MEX files (to call
custom C or Fortran routines directly from Matlab) can be executed in Octave
with some limitations.


Octave runs are typically slower than the equivalent Matlab run. However,
Octave is not license-limited, so many more simultaneous runs are possible than
with Matlab. Thus, on the Biowulf cluster Octave is most useful for projects
which can be split into large numbers of independent simultaneous runs.



Documentation
* Online documentation is available by typing help at the Octave
prompt. (e.g. help acos)
* [The Octave
Manual](http://www.gnu.org/software/octave/doc/interpreter/)
* [Octave Wiki](http://wiki.octave.org/)
* [Differences between Matlab & Octave](http://en.wikibooks.org/wiki/MATLAB_Programming/Differences_between_Octave_and_MATLAB), at Wikibooks.
* [Efficient Matlab and Octave](http://homepages.inf.ed.ac.uk/imurray2/compnotes/matlab_octave_efficiency.html) -- by Iain Murray, University of Edinburgh.


Important Notes
* Module Name: octave (see [the modules page](/apps/modules.html) for more information)
* Multithreaded/singlethreaded/MPI...
* The following libraries are currently not available with Octave on Biowulf. 


|  |  |
| --- | --- |
| OSMesa library:  | Offscreen rendering with OpenGL will be disabled. |
| Qscintilla library:  |  Built-in GUI editor is disabled. |



Interactive session with GUI
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


[user@cn3144 ~]$ **module load octave**
[+] Loading gnuplot 5.2.2  ...
[+] Loading HDF5  1.10.1
[+] Loading LAPACK/3.8.0  3.8.0-gcc4.8.5  libraries...
[+] Loading FFTW 3.3.7 , compiled with gcc4.8.5  and openmpi2.1.2  ...
[+] Loading ImageMagick  7.0.7  on cn4242
[+] Loading gl2ps 1.4.0  ...
[+] Loading Octave 4.4.0  ...

[user@cn3144 ~]$ **octave --gui**
![](/images/octave1.jpg)

```


Interactive session without GUI
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

[user@cn3144 ~]$ **module load octave**
[+] Loading gnuplot 5.2.2  ...
[+] Loading HDF5  1.10.1
[+] Loading LAPACK/3.8.0  3.8.0-gcc4.8.5  libraries...
[+] Loading FFTW 3.3.7 , compiled with gcc4.8.5  and openmpi2.1.2  ...
[+] Loading ImageMagick  7.0.7  on cn4242
[+] Loading gl2ps 1.4.0  ...
[+] Loading Octave 4.4.0  ...

[user@cn3144 ~]$ **octave-cli**
warning: function ./test.m shadows a core library function
GNU Octave, version 4.0.2
Copyright (C) 2016 John W. Eaton and others.
This is free software; see the source code for copying conditions.
There is ABSOLUTELY NO WARRANTY; not even for MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  For details, type 'warranty'.

Octave was configured for "x86_64-unknown-linux-gnu".

Additional information about Octave is available at http://www.octave.org.

Please contribute if you find this software useful.
For more information, visit http://www.octave.org/get-involved.html

Read http://www.octave.org/bugs.html to learn how to submit bug reports.
For information about changes from previous versions, type 'news'.

octave:1>**a=1+1**
a =  2
octave:2>  **sin(a)**
ans = -0.27942
octave:3> **a=5\*pi**
x =  15.708
octave:4> **sin(a)**
ans =  6.1230e-16
octave:6> **quit**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).

1. Create an Octave script. The following sample script takes a single
argument from the command line.




```

!/bin/env octave -qf

if( nargin != 1 )
    printf( "Usage: %s <integer>\n", program_name );
    return;
endif

len = str2num( argv(){1} );
printf( "Working with array size %6d\n", len );

clear a; 
tic(); 
for i=1:len
    a(i) = i; 
endfor 
time1 = toc();

a = [1]; 
tic(); 
for i=2:len 
    a = [a i]; 
endfor
time2 = toc();

printf( "The time taken for method 1 was %.4f seconds\n", time1 );
printf( "The time taken for method 2 was %.4f seconds\n", time2 );

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.


- Create a batch script for the job. Sample script:


```

-------- file array.bat----------------------
#!/bin/bash

module load octave/4.0.2

cd /data/user/mydir
./array.oct 10000 
---------------------------------------------

```

- Submit this job to the batch system with:

```

sbatch  [--mem=#g] array.bat

```


Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. octave.swarm). For example:



```

array.oct 1000 > array.out.1000
array.oct 2000 > array.out.2000
array.oct 3000 > array.out.3000
[...]

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f octave.swarm [-g #] [-t #] --module octave
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module octave Loads the octave module for each subjob in the swarm 
 | |
 | |
 | |



Octave Packages
If you have a particular package from the [Octave-Forge project](http://octave.sourceforge.net/packages.php) that you need installed, let us know. To see the centrally installed packages,
use pkg list e.g.

```

octave:3> **pkg list**
Package Name  | Version | Installation directory
--------------+---------+-----------------------
     control  |   2.8.2 | /usr/local/apps/octave/4.0.2/share/octave/packages/control-2.8.2
       image  |   2.0.0 | /usr/local/apps/octave/4.0.2/share/octave/packages/image-2.0.0
      signal  |   1.3.2 | /usr/local/apps/octave/4.0.2/share/octave/packages/signal-1.3.2


```


You can also install your
own Octave tool-boxes in your home directory with pkg install. You will need to set the TMPDIR variable first. e.g.

```


cn1234% **export TMPDIR=/lscratch/$USER**

cn1234% **octave**
octave:1> **pkg install -forge fpl**
For information about changes from previous versions of the fpl package, run 'news fpl'.

octave:2> **pkg list**
Package Name  | Version | Installation directory
--------------+---------+-----------------------
     control  |   2.8.2 | /usr/local/apps/octave/4.0.2/share/octave/packages/control-2.8.2
         fpl  |   1.3.5 | /home/teacher/octave/fpl-1.3.5
       image  |   2.0.0 | /usr/local/apps/octave/4.0.2/share/octave/packages/image-2.0.0
      signal  |   1.3.2 | /usr/local/apps/octave/4.0.2/share/octave/packages/signal-1.3.2

```


















