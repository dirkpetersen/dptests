

document.querySelector('title').textContent = ' Neuron on Biowulf';

Neuron on Biowulf



|  |
| --- |
| 
Quick Links
[Batch job on Biowulf](#batch)
[Interactive job on Biowulf](#int)
[Python Interface](#python)
[Documentation](#doc)
 |



[![](/images/neuron.jpg)](http://www.neuron.yale.edu/neuron/)

NEURON is a simulation environment for modeling individual neurons and networks of neurons. It provides tools for conveniently building, managing, and using models in a way that is numerically sound and computationally efficient. It is particularly well-suited to problems that are closely linked to experimental data, especially those that involve cells with complex anatomical and biophysical properties.

Neuron was developed at Yale University. [[Neuron website](http://www.neuron.yale.edu)]

Neuron can be run via the GUI, which is best done on Helix. The purpose of running Neuron on Biowulf is to run batch jobs using command-line scripts, which can be multithreaded or parallel.


Documentation

* [Neuron docs on the Yale site](http://www.neuron.yale.edu/neuron/docs).
* [Neuron courses](http://www.neuron.yale.edu/neuron/courses).




interactive job

As a first test, try 'neurondemo'. After some preparatory command-line output, you should see several Neuron windows appear on your screen. 

```

[user@biowulf]$ **sinteractive**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **neurondemo**

![](/images/neuron1.jpg)
![](/images/neuron2.jpg) ![](/images/neuron3.jpg)![](/images/neuron4.jpg)

[user@cn3144 ~]$ **exit**






```


To make Neuron read your input file 'foo.hoc', type 'nrngui foo.hoc'.

Neuron can potentially multithread -- use multiple cores on Helix. However, this will only help if your simulations have cells with lots of compartments and mechanisms, i.e. requiring at least a couple thousand ODEs per cell? If no, multithreaded won't help. It is recommended that you start by using the GUI RunControl panel and the Neuron Main Menu -> Tools -> Parallel Computing tool -- to run a few interactive tests to decide on appropriate simulation parameters. See [for more info](http://www.neuron.yale.edu/phpbb/viewtopic.php?f=31&t=1710).



Batch job on Biowulf

Create a batch script along the following lines:

```

#/bin/bash

module load neuron/7.5-mpi
mpirun nrniv -mpi myfile.hoc

```


where myfile.hoc is the file containing the hoc source code. 

Submit this job with:

```

biowulf% sbatch --ntasks=# myfile.bat

```

This will submit the job to the number of cores specified by the '--ntasks=#' parameter. OpenMPI will obtain the number of processes and the list of nodes/cores from Slurm directly, so you don't have to specify those parameters in the mpirun command in your batch script. If you need more than the default 1 GB of memory per core, you should specify that with

```

biowulf% sbatch --ntasks=# --mem-per-cpu=4g myfile.bat

```



Python interface to Neuron

The Python language interface for Neuron is installed in v 7.4 (note: not in the version 7.4-mpi at present). Sample session:


```

[user@biowulf ~]$ **sinteractive**
salloc.exe: Pending job allocation 18431473
salloc.exe: job 18431473 queued and waiting for resources
salloc.exe: job 18431473 has been allocated resources
salloc.exe: Granted job allocation 18431473
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn1755 are ready for job

[user@cn1755 ~]$ **module load neuron/7.4**
[+] Loading python 2.7.6 ...
[+] Loading Neuron 7.4
 ...
[user@cn1755 ~]$ **ipython**
Python 2.7.6 |Continuum Analytics, Inc.| (default, May 27 2014, 14:50:58)
Type "copyright", "credits" or "license" for more information.

IPython 4.2.0 -- An enhanced Interactive Python.
?         -> Introduction and overview of IPython's features.
%quickref -> Quick reference.
help      -> Python's own help system.
object?   -> Details about 'object', use 'object??' for extra details.

In [1]: import neuron
NEURON -- Release 7.4 (1370:16a7055d4a86) 2015-11-09
Duke, Yale, and the BlueBrain Project -- Copyright 1984-2015
See http://www.neuron.yale.edu/neuron/credits

[.....]

In [20]: **exit**
[user@cn1755 ~]$ exit

```







































