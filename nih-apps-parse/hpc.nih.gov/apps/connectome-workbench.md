

document.querySelector('title').textContent = ' Connectome Workbench on Biowulf';

Connectome Workbench on Biowulf



|  |
| --- |
| 
Quick Links
[Interactive job on Biowulf](#int)
[Batch job on Biowulf](#serial)
[Swarm of jobs](#swarm)
[Documentation](#doc)
 |



Description
Connectome Workbench is an open source, freely available visualization and 
discovery tool used to map neuroimaging data, especially data generated by the 
Human Connectome Project. 


### Reference


* Van Essen, D. C., Smith, S. M., Barch, D. M., Behrens, T. E., 
 Yacoub, E., Ugurbil, K., & WU-Minn HCP Consortium. (2013). The WU-Minn human 
 connectome project: an overview. *Neuroimage*, *80*, 62-79.
* Glasser, M. F., Sotiropoulos, S. N., Wilson, J. A., Coalson, T. S., 
 Fischl, B., Andersson, J. L., ... & WU-Minn HCP Consortium. (2013). The 
 minimal preprocessing pipelines for the Human Connectome Project. 
 *Neuroimage*, *80*, 105-124.
* Marcus, D. S., Harwell, J., Olsen, T., Hodge, M., Glasser, M. F., Prior, 
 F., ... & Van Essen, D. C. (2011). Informatics and data mining tools and 
 strategies for the human connectome project. *Frontiers in 
 neuroinformatics*, *5*.


### Web sites


* [Home page](http://www.humanconnectome.org/software/connectome-workbench.html)
* [GitHub](https://github.com/Washington-University/workbench)



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

[user@cn3144 ~]$ **module load connectome-workbench**

[user@cn3144 ~]$ **wb\_command -help**

[user@cn3144 ~]$ **wb\_view -help**

[user@cn3144 ~]$ **wb\_import**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Additionally, if you have X11 forwarding enabled you can use the workbench gui
by entering **wb\_view** with no arguments.





Batch jobs on Biowulf

wb\_command gives access to a wide variety of commands related to data analysis 
iand image processing that can be used in your scripts. Use the following 
argument to view the documentation.




```

[user@helix ~]$ **wb\_command -all-commands-help**

```


After selecting the command(s) you want to run in a batch script, set it up like 
so.




```

#!/bin/bash
# this file is called myjob.bat

cd /data/$USER/my_workbench_data
module load connectome-workbench

wb_command ...my_command1...
wb_command ...my_command2...
wb_command ...my_command3...
...
wb_command ...my_commandN...

```


Submit the job to SLURM with:




```

[user@helix ~]$ **sbatch myjob.bat**

```



Swarm of jobs on Biowulf

If you need to run a large number of batch jobs, you can run them in parallel 
using the swarm program. For this example, imagine we have a large number of 
files to convert from cifti to nifti format. We would create a swarm file 
(called "workbench.swarm") that looked something like this:




```

 wb_command -cifti-convert -to-nifti /path/to/file1.cifti /path/to/file1.nifti
 wb_command -cifti-convert -to-nifti /path/to/file2.cifti /path/to/file2.nifti
 wb_command -cifti-convert -to-nifti /path/to/file3.cifti /path/to/file3.nifti
 ...                             
 wb_command -cifti-convert -to-nifti /path/to/fileN.cifti /path/to/fileN.nifti

```

We would submit it to the queue with



```

[user@helix ~]$ **swarm -f workbench.swarm --module connectome-workbench** 

```


See the [swarm webpage](/apps/swarm.html) for more
information, or contact the Biowulf staff at staff@hpc.nih.gov 






Documentation
* [Connectome Workbench Tutorial](http://www.humanconnectome.org/documentation/tutorials/Connectome_WB_Tutorial_v1.0.pdf)
* [Using Workbench Command](http://www.humanconnectome.org/software/workbench-command.php)





