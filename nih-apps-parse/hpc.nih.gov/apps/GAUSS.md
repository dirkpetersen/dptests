

document.querySelector('title').textContent = 'GAUSS on Biowulf';
GAUSS on Biowulf


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



![gauss](/images/gauss.jpg)**[The GAUSS Mathematical
and Statistical System](http://www.aptech.com/products/gauss-mathematical-and-statistical-system/)** is a fast matrix programming language widely
used by scientists, engineers, statisticians, biometricians, econometricians,
and financial analysts. Designed for computationally intensive tasks, the GAUSS
system is ideally suited for the researcher who does not have the time required
to develop programs in C or FORTRAN but finds that most statistical or
mathematical "packages" are not flexible or powerful enough to perform
complicated analysis or to work on large problems.
The GAUSS executables are not multithreaded or parallel. The advantage of
using GAUSS on Biowulf would be to run many GAUSS jobs simultaneously, i.e. a
'swarm' of single-threaded jobs.






Documentation
* GAUSS has online help. Typing 'gauss' at the Biowulf prompt will bring up
the GAUSS Xwindows interface. Click on the Help button to view the GAUSS
help.
* [GAUSS. an introduction](http://personal.lse.ac.uk/goujard/teaching0809/EC475/gauss_introduction.pdf), at the London School of Economics


Important Notes
* Module Name: gauss (see [the modules page](/apps/modules.html) for more information)
* There are 2 versions of GAUSS available, as shown in the table below. There are a limited
number of licenses for GAUSS 10, so you will need to use GAUSS 3.2 to run swarm jobs.


|  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| **Version **To use **Licenses **Architecture **Can be run on
| Gauss 10 module load gauss/10 2   64-bit  any Biowulf node
| Gauss 3.2 module load gauss/3.2 unlimited  32-bit  any Biowulf node
 | | | | |
 | | | | |** |** |** |** |** |



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

[user@cn3144 ~]$ cat > doloop.inp
format /rdn 1,0; 
space = "      "; 
comma = ","; 
i = 1;
do while i <= 4; 
j = 1;
do while j <= 3; print space i comma j;; 
j = j+1;
endo; 
i = i+1; 
print;
endo;

(type Ctrl-D)

[user@cn3144 ~]$ **module load gauss**
[+] Loading GAUSS  10  on cn3144

[user@cn3144 ~]$ **tgauss doloop.inp**

GAUSS 10.0.3 (Dec 22 2009, 1346) 64-bit
(C)Copyright 1984-2009 Aptech Systems, Inc.
All Rights Reserved Worldwide.

      1,1      1,2      1,3
      2,1      2,2      2,3
      3,1      3,2      3,3
      4,1      4,2      4,3
      
(gauss) quit

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. GAUSS.sh). For example:



```

#!/bin/tcsh

module load gauss/3.2
cd mydir
gauss -v -b gauss.in > gauss.out


```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch  [--mem=#] GAUSS.sh
```


### Multi-threaded GAUSS jobs



There are several threading programs that can be used to multithread (i.e. run on multiple processors) specific parts of your programs. These are described [here](http://www.aptech.com/resources/tutorials/) . These threading functions can be used to utilize all the processors on an allocated node. It is important to know exactly how many threads you are executing and match this number to the available processors on the node, so that you neither overload the node (very inefficient) or waste processors. 

For example, the following sample code from the GAUSS 10 User Guide defines 4 concurrent threads:

```

ThreadStat n = m'm;    //Thread 1
ThreadBegin;           //Thread 2  
y = x'x;
z = y'y; 
ThreadEnd;
ThreadBegin;          //Thread 3
q = r'r; 
r = q'q;
ThreadEnd; 
ThreadStat p = o'o;   //Thread 4

```


Write a batch script along the following lines:

```

#!/bin/bash
# ---- this file is called myjobscript -------

module load gauss/10
cd mydir
tgauss multi.inp

```


This job can be submitted with:

```

sbatch --cpus-per-task=4  [--mem=#g] myjobscript

```

This command will submit the job to 2 cores (4 CPUS to match the 4 threads) [and # GB of memory if specified]. 


Swarm of GAUSS 3.2 Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.

Note: there are unlimited GAUSS 3.2 licenses, so you can only use v3.2 to run a swarm of jobs. 
Create a swarmfile (e.g. GAUSS.swarm). For example:



```

gauss -v -b gauss1.in > gauss1.out
gauss -v -b gauss2.in > gauss2.out
gauss -v -b gauss3.in > gauss3.out
gauss -v -b gauss4.in > gauss4.out
gauss -v -b gauss5.in > gauss5.out

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f GAUSS.swarm [-g #] [-t #] --module gauss/3.2
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module GAUSS Loads the GAUSS module for each subjob in the swarm 
 | |
 | |
 | |






















