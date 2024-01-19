

document.querySelector('title').textContent = 'autoreject on Biowulf';
autoreject on Biowulf


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



This is a library to automatically reject bad trials and repair bad sensors in magneto-/electroencephalography (M/EEG) data.




![autoreject](/images/autoreject.png)



### References:


* Mainak Jas et.al.*Autoreject: Automated artifact rejection for MEG and EEG data*
NeuroImage, 2017 159, 417-429.
 [PubMed](https://www.ncbi.nlm.nih.gov/pubmed/28645840) | 
 [Journal](https://www.sciencedirect.com/science/article/pii/S1053811917305013)


Documentation
* autoreject Main Site:[Main Site](https://autoreject.github.io/stable/index.html)


Important Notes
* Module Name: autoreject (see [the modules page](/apps/modules.html) for more information)
	autoreject is installed as a singularity container, which means it came with it's own python environment. 
	
	
	It could be run as command line or through jupyter notebook.* Environment variables set 
	+ $AUTOREJECT\_TEST\_DATA



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive --cpus-per-task=10 --mem=10G --gres=lscratch:200 --tunnel**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job
Created 1 generic SSH tunnel(s) from this compute node to                  
biowulf for your use at port numbers defined                               
in the $PORTn ($PORT1, ...) environment variables.                         
                                                                           
                                                                           
Please create a SSH tunnel from your workstation to these ports on biowulf.
On Linux/MacOS, open a terminal and run:                                   
                                                                           
    ssh  -L 33327:localhost:33327 biowulf.nih.gov                          
                                                                           
For Windows instructions, see https://hpc.nih.gov/docs/tunneling          
[user@cn3144]$ **module load autoreject**
[user@cn3144]$ **python -c "import autoreject"**
[user@cn3144]$ **cd /lscratch/$SLURM\_JOB\_ID**
[user@cn3144]$ **cp -r ${AUTOREJECT\_TEST\_DATA:-none}/mne\_data .**
[user@cn3144]$ **jupyter notebook --ip localhost --port $PORT1 --no-browser** 
[I 17:11:40.505 NotebookApp] Serving notebooks from local directory
[I 17:11:40.505 NotebookApp] Jupyter Notebook 6.4.10 is running at:
[I 17:11:40.505 NotebookApp] http://localhost:37859/?token=xxxxxxxx
[I 17:11:40.506 NotebookApp]  or http://127.0.0.1:37859/?token=xxxxxxx
[I 17:11:40.506 NotebookApp] Use Control-C to stop this server and shut down all kernels (twice to skip confirmation).
[C 17:11:40.512 NotebookApp]

    To access the notebook, open this file in a browser:
        file:///home/apptest1/.local/share/jupyter/runtime/nbserver-29841-open.html
    Or copy and paste one of these URLs:
        http://localhost:37859/?token=xxxxxxx
     or http://127.0.0.1:37859/?token=xxxxxxx

```

 Then you can open a browser from your computer to connect to the jupyter notebook:
 ![autoreject_jupyter](/images/autoreject_jupyter.png)



```

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. autoreject.sh). For example:



```

#!/bin/bash
#SBATCH --cpus-per-task=10
#SBATCH --mem=10G
#SBATCH --time=2:00:00
#SBATCH --partition=norm

set -e
module load autoreject
# command 

```

 Submit the job:

```
sbatch autoreject.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. job.swarm). For example:



```

cmd1
cmd2
cmd3

    
```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f job.swarm [-g #] --module autoreject
```

where
 

|  |  |
| --- | --- |
| -g *#*  | Number of Gigabytes of memory required for each process
 (1 line in the swarm command file)  |
| --module  | Loads the module for each subjob in the swarm  |












