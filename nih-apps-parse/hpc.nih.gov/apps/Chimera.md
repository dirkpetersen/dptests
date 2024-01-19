

document.querySelector('title').textContent = 'Chimera on Biowulf';
Chimera on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
 |


UCSF Chimera is a highly extensible program for interactive visualization and analysis of molecular structures and related data, including density maps, supramolecular assemblies, sequence alignments, docking results, trajectories, and conformational ensembles. High-quality images and animations can be generated. Chimera includes complete documentation and several tutorials, and can be downloaded free of charge for academic, government, non-profit, and personal use.



### References:


* Please see <https://www.cgl.ucsf.edu/chimera/docs/credits.html> for details on how to cite UCSF Chimera.


Documentation
* [Chimera Documentation](http://www.cgl.ucsf.edu/chimera/docindex.html)
* [Chimera Programmer's Guide](http://www.cgl.ucsf.edu/chimera/docs/ProgrammersGuide/index.html)


Important Notes
* Module Name: Chimera (see [the modules page](/apps/modules.html) for more information)
* Multithreaded/singlethreaded/MPI...
* environment variables set 
	+ LD\_LIBRARY\_PATH
	+ PYTHONPATH* Reference data in /pdb



This application requires an [X-Windows connection](/docs/connect.html). 


There are combinations of X11 servers and drivers that cause Chimera to crash. It is known that XQuartz (v2.7.x) is incompatible with Chimera. Users are encouraged to use [NX or FastX](https://hpc.nih.gov/docs/nx.html) as their X11 servers.


Chimera will not run normally on GPU nodes. Chimera can only run on CPU-only nodes


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

[user@cn3144 ~]$ module load Chimera
[user@cn3144 ~]$ chimera &

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```

Running a Chimera script

A Chimera script contains a series of commands for Chimera. Here is an example of how to run a demo script from the [Chimera Scripting Tutorial](https://www.cgl.ucsf.edu/Outreach/Workshops/UCSF-Fall-2005/09-ScriptDemo/Scripting_Tutorial.html).

Create a file with the Chimera commands you want to run. In this case, the demo commands will load the Green Fluorescent Protein and display it in various ways.
This sample script is available in /usr/local/apps/Chimera/examples/GFP.chimera. 

Open a [NX session to Biowulf](/docs/nx.html).


```

biowulf% **sinteractive** 					#Start an interactive sesssion on Biowulf
salloc: Pending job allocation 48280698
salloc: job 48280698 queued and waiting for resources
salloc: job 48280698 has been allocated resources
salloc: Granted job allocation 48280698
salloc: Waiting for resource configuration
salloc: Nodes cn4271 are ready for job

[user@cn4271] **cp /usr/local/apps/Chimera/examples/GFP.chimera /home/$USER**

[user@cn4271]  **module load chimera**

[user@cn4271]  **chimera GFP.chimera &**




Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).

Chimera is almost entirely written in python. As such, scripting chimera functions is relatively straightforward. Once a python script is written, it can be run by either including this shebang at the top, 



 
```
#!/usr/bin/env python2.7
```


 changing the script to executable,



 
```
chmod +x myChimeraScript.py
```


 and then running it like so:



 
```
./myChimeraScript.py
```


 Or, simply leaving off the shebang and calling the correct python executable:



 
```
python2.7 myChimeraScript.py
```


 [Here](http://www.cgl.ucsf.edu/chimera/docs/ProgrammersGuide/Examples/CreateMolecule.py) is an 
 example of a chimera python script.



Create a batch input file (e.g. Chimera.sh).




```

#!/bin/bash
module load Chimera
./myChimeraScript.py

```


Submit this job using the Slurm [sbatch](/docs/userguide.html) command.




```
sbatch [--cpus-per-task=#] [--mem=#] Chimera.sh
```








```










