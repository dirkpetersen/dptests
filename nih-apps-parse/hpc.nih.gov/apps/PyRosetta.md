

document.querySelector('title').textContent = 'PyRosetta on Biowulf';
PyRosetta on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
 |



PyRosetta is an interactive Python-based interface to the powerful Rosetta molecular modeling suite. It enables users to design their own custom molecular modeling algorithms using Rosetta sampling methods and energy functions.




PyRosetta was created at Johns Hopkins University by Jeffrey J. Gray, Sergey Lyskov, and the PyRosetta Team.



Documentation
* [PyRosetta Main Site](http://www.pyrosetta.org/)
* [PyRosetta Tutorial Workshop](http://www.pyrosetta.org/tutorials)


Important Notes
* Module Name: PyRosetta (see [the modules page](/apps/modules.html) for more information)
* Singlethreaded
* Environment variables set 
	+ PYROSETTA\_HOME* Example files in $PYROSETTA\_HOME/demo* Reference data in /pdb/pdb


Many aspects of PyRosetta require an X11 connection. Please see <https://hpc.nih.gov/docs/connect.html> for graphical methods of connecting.


Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive --gres=lscratch:10**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load PyRosetta**
[user@cn3144 ~]$ **python**
Python 2.7.15 |Anaconda, Inc.| (default, Oct 23 2018, 18:31:10) 
[GCC 7.3.0] on linux2
Type "help", "copyright", "credits" or "license" for more information.
>>> **import pyrosetta**
>>> **pyrosetta.init()**

... your code here ...

**quit()**
[user@cn3144 ~]$ **cd /lscratch/$SLURM\_JOB\_ID**
[user@cn3144 ~]$ **cp -r $PYROSETTA\_HOME/demo .**
[user@cn3144 ~]$ **cd demo**
[user@cn3144 ~]$ **cp -r $PYROSETTA\_HOME/test .**
[user@cn3144 ~]$ **mkdir .test.output**
[user@cn3144 ~]$ **python D010\_Pose\_structure.py**

... pyrosetta output here ...

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. PyRosetta.sh). For example:



```

#!/bin/bash
module load PyRosetta
python your_python_code.py

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] PyRosetta.sh
```





