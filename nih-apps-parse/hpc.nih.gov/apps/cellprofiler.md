

document.querySelector('title').textContent = 'cellprofiler on Biowulf';
cellprofiler on Biowulf


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



CellProfiler is a free open-source software designed to enable biologists without training in computer vision or programming to quantitatively measure phenotypes from thousands of images automatically.



Documentation
* [cellprofiler Main Site](https://github.com/CellProfiler/CellProfiler)


Important Notes
* Module Name: cellprofiler (see [the modules page](/apps/modules.html) for more information)
* singlethreaded (non-GUI/headless version)
* environment variables set 
	+ CELLPROFILER\_TESTDATA=/usr/local/apps/cellprofiler/TESTDATA* Example files in /usr/local/apps/cellprofiler/TESTDATA



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf ~]$ **sinteractive -c2 --mem=4g --gres=lscratch:10**
salloc.exe: Pending job allocation 11258266
salloc.exe: job 11258266 queued and waiting for resources
salloc.exe: job 11258266 has been allocated resources
salloc.exe: Granted job allocation 11258266
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn0884 are ready for job
srun: error: x11: no local DISPLAY defined, skipping
error: unable to open file /tmp/slurm-spank-x11.11258266.0
slurmstepd: error: x11: unable to read DISPLAY value

[user@cn0884 ~]$ **cd /lscratch/$SLURM\_JOB\_ID**

[user@cn0884 11258266]$ **module load cellprofiler**
[+] Loading cellprofiler 3.1.8  ...

[user@cn0884 11258266]$ **cp ${CELLPROFILER\_TESTDATA}/\* .**

[user@cn0884 11258266]$ **unzip ExampleHuman.zip**
Archive:  ExampleHuman.zip
   creating: ExampleHuman/
  inflating: ExampleHuman/ExampleHuman.cppipe
   creating: __MACOSX/
   creating: __MACOSX/ExampleHuman/
  inflating: __MACOSX/ExampleHuman/._ExampleHuman.cppipe
   creating: ExampleHuman/images/
  inflating: ExampleHuman/images/AS_09125_050116030001_D03f00d0.tif
   creating: __MACOSX/ExampleHuman/images/
  inflating: __MACOSX/ExampleHuman/images/._AS_09125_050116030001_D03f00d0.tif
  inflating: ExampleHuman/images/AS_09125_050116030001_D03f00d1.tif
  inflating: __MACOSX/ExampleHuman/images/._AS_09125_050116030001_D03f00d1.tif
  inflating: ExampleHuman/images/AS_09125_050116030001_D03f00d2.tif
  inflating: __MACOSX/ExampleHuman/images/._AS_09125_050116030001_D03f00d2.tif
  inflating: __MACOSX/ExampleHuman/._images
  inflating: ExampleHuman/README.md
  inflating: __MACOSX/ExampleHuman/._README.md
  inflating: __MACOSX/._ExampleHuman

[user@cn0884 11258266]$ **cd ExampleHuman/**

[user@cn0884 ExampleHuman]$ **cellprofiler -p ExampleHuman.cppipe -c -r -i images/**
/usr/local/Anaconda/envs_app/cellprofiler/3.1.8/lib/python2.7/site-packages/cellprofiler/utilities/hdf5_dict.py:539: FutureWarning: Conversion of the second argument of issubdtype from `int` to `np.signedinteger` is deprecated. In future, it will be treated as `np.int64 == np.dtype(int).type`.
  np.issubdtype(hdf5_type, int) or
/usr/local/Anaconda/envs_app/cellprofiler/3.1.8/lib/python2.7/site-packages/cellprofiler/utilities/hdf5_dict.py:541: FutureWarning: Conversion of the second argument of issubdtype from `float` to `np.floating` is deprecated. In future, it will be treated as `np.float64 == np.dtype(float).type`.
  hdf5_type_is_float = np.issubdtype(hdf5_type, float)
Times reported are CPU and Wall-clock times for each module
Wed Mar 24 13:33:33 2021: Image # 1, module Images # 1: CPU_time = 0.00 secs, Wall_time = 0.00 secs
Wed Mar 24 13:33:33 2021: Image # 1, module Metadata # 2: CPU_time = 0.00 secs, Wall_time = 0.00 secs
Wed Mar 24 13:33:33 2021: Image # 1, module NamesAndTypes # 3: CPU_time = 1.87 secs, Wall_time = 0.97 secs
Wed Mar 24 13:33:34 2021: Image # 1, module Groups # 4: CPU_time = 0.00 secs, Wall_time = 0.00 secs
/usr/local/Anaconda/envs_app/cellprofiler/3.1.8/lib/python2.7/site-packages/centrosome/cpmorphology.py:4209: FutureWarning: Using a non-tuple sequence for multidimensional indexing is deprecated; use `arr[tuple(seq)]` instead of `arr[seq]`. In the future this will be interpreted as an array index, `arr[np.array(seq)]`, which will result either in an error or a different result.
  big_labels[[slice(fe,-fe) for fe in footprint_extent]] = labels
[...snip]
Wed Mar 24 13:33:41 2021: Image # 1, module SaveImages # 13: CPU_time = 0.35 secs, Wall_time = 0.20 secs
Wed Mar 24 13:33:41 2021: Image # 1, module ExportToSpreadsheet # 14: CPU_time = 0.00 secs, Wall_time = 0.00 secs

[user@cn0884 ExampleHuman]$ **exit**
exit
salloc.exe: Relinquishing job allocation 11258266
salloc.exe: Job allocation 11258266 has been revoked.

[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. cellprofiler.sh). For example:



```

#!/bin/bash
set -e
module load cellprofiler
cellprofiler -p ExampleHuman.cppipe -c -r -i images/

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] cellprofiler.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. cellprofiler.swarm). For example:



```

cellprofiler -p ExampleHuman1.cppipe -c -r -i images1/
cellprofiler -p ExampleHuman2.cppipe -c -r -i images2/
cellprofiler -p ExampleHuman3.cppipe -c -r -i images3/

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f cellprofiler.swarm [-g #] [-t #] --module cellprofiler
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module cellprofiler Loads the cellprofiler module for each subjob in the swarm 
 | |
 | |
 | |








