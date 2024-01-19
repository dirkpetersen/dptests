

document.querySelector('title').textContent = 'AMPL on Biowulf';
AMPL on Biowulf


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



AMPL is the Accelerating Therapeutics for Opportunites in Medicine (ATOM) Consortium Modeling PipeLine for drug discovery.



### Reference:


* [Minnich, Amanda J., et al. "AMPL: A Data-Driven Modeling Pipeline for Drug Discovery." *Journal of Chemical Information and Modeling* 60.4 (2020): 1955-1968.](https://pubs.acs.org/doi/abs/10.1021/acs.jcim.9b01053)


Documentation
* [AMPL Read the Docs site](https://ampl.readthedocs.io/en/latest/pipeline.html)
* [AMPL Tutorial](https://github.com/ATOMconsortium/AMPL/tree/master/atomsci/ddm/examples/tutorials)


Important Notes
* Module Name: ampl (see [the modules page](/apps/modules.html) for more information)
 * The ampl module sets the $AMPL\_HOME environment variable which points to the location of the installation. A matched version of the GitHub repo is present in this directory for copying. If you want to modify ampl scripts before running, you should copy this directory to your own space. 
 * AMPL is installed within a pre-made conda environment. All binaries in the original environment are symlinked into your path (it is mostly overwriting your existing Python). Do not attempt to use simultaneously with another Python environment.
 * If you find that your sys.path contains libraries within your home directory that conflict with the central installation of ampl, you can start python with the -Es flags, or execute export PYTHONNOUSERSITE=True.
 * The following commands and scripts are avaiable after loading the module. Please note that some of these scripts may not be executable. Please contact the NIH HPC staff if other commands are required within the AMPL environment:
	+ ave\_splitter.py
	+ chem\_diversity.py
	+ compare\_models.py
	+ curate\_data.py
	+ data\_curation\_functions.py
	+ datastore\_functions.py
	+ dist\_metrics.py
	+ diversity\_plots.py
	+ featurization.py
	+ genTestset.py
	+ hyperparam\_search\_wrapper.pec/ampl.sh
	+ hyper\_perf\_plots.py
	+ ipython
	+ jupyter
	+ jupyter-notebook
	+ jupyter-run
	+ llnl\_utils.py
	+ model\_datasets.py
	+ model\_pipeline.py
	+ model\_tracker.py
	+ model\_wrapper.py
	+ nosetests
	+ open-docs.py
	+ parameter\_parser.py
	+ perf\_data.py
	+ perf\_plots.py
	+ process\_slurm.py
	+ pubchem\_utils.py
	+ pytest
	+ python
	+ rdkit\_easy.py
	+ splitting.py
	+ struct\_utils.py
	+ temporal\_splitter.py
	+ tensorboard
	+ transformations.py



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf ~]$ **sinteractive --gres=lscratch:10 --mem=8g -c4 --tunnel**
salloc: Pending job allocation 35296239
salloc: job 35296239 queued and waiting for resources
salloc: job 35296239 has been allocated resources
salloc: Granted job allocation 35296239
salloc: Waiting for resource configuration
salloc: Nodes cn0857 are ready for job
srun: error: x11: no local DISPLAY defined, skipping
error: unable to open file /tmp/slurm-spank-x11.35296239.0
slurmstepd: error: x11: unable to read DISPLAY value

Created 1 generic SSH tunnel(s) from this compute node to
biowulf for your use at port numbers defined
in the $PORTn ($PORT1, ...) environment variables.


Please create a SSH tunnel from your workstation to these ports on biowulf.
On Linux/MacOS, open a terminal and run:

    ssh  -L 36858:localhost:36858 user@biowulf.nih.gov

For Windows instructions, see https://hpc.nih.gov/docs/tunneling


[user@cn0857 ~]$ **module load ampl**
[+] Loading ampl  1.3.0  on cn0857

[user@cn0857 ~]$ **ipython**
Python 3.7.10 (default, Feb 26 2021, 18:47:35)
Type 'copyright', 'credits' or 'license' for more information
IPython 7.16.1 -- An enhanced Interactive Python. Type '?' for help.

In [1]: **import atomsci**

In [2]: **quit()**

[user@cn0857 ~]$ **cd /lscratch/$SLURM\_JOB\_ID**

[user@cn0857 35296239]$ **cp -r ${AMPL\_HOME}/ampl .**

[user@cn0857 35296239]$ **cd ampl/atomsci/ddm/**

[user@cn0857 ddm]$ **jupyter-notebook --no-browser --port=$PORT1**
[I 14:23:09.610 NotebookApp] Serving notebooks from local directory: /lscratch/35296239/ampl/atomsci/ddm
[I 14:23:09.610 NotebookApp] Jupyter Notebook 6.4.10 is running at:
[I 14:23:09.610 NotebookApp] http://localhost:36858/?token=5ebfcff3cfccd637d36250f9a5a8191a22011edd2a262db5
[I 14:23:09.610 NotebookApp]  or http://127.0.0.1:36858/?token=5ebfcff3cfccd637d36250f9a5a8191a22011edd2a262db5
[I 14:23:09.610 NotebookApp] Use Control-C to stop this server and shut down all kernels (twice to skip confirmation).
[C 14:23:09.618 NotebookApp]

    To access the notebook, open this file in a browser:
        file:///spin1/home/linux/user/.local/share/jupyter/runtime/nbserver-63375-open.html
    Or copy and paste one of these URLs:
        http://localhost:36858/?token=5ebfcff3cfccd637d36250f9a5a8191a22011edd2a262db5
     or http://127.0.0.1:36858/?token=5ebfcff3cfccd637d36250f9a5a8191a22011edd2a262db5

```

Now in a new terminal on your local machine, you should initiate a new ssh session to Biowulf to establish the ssh tunnel. Using the port that was allocated to you during the sinteractive step above:

```

[user@mymachine ~]$ **ssh -L 36858:localhost:36858 user@biowulf.nih.gov**

[user@biowulf ~]$


```

This will allow you to open a local browser and view the Jupyter Notebook session running on the compute node:

 ![AMPL Jupyter Notebook image](/images/ampl-jupyter-example.PNG)



When you are finished close the Jupyter Notebook session in original terminal using Ctrl-C and exit the interactive session to conserve system resources. 

```

**^C**[I 11:48:01.862 NotebookApp] interrupted
Serving notebooks from local directory: /lscratch/1151805/ampl/atomsci/ddm
1 active kernel
The Jupyter Notebook is running at:
http://localhost:36858/?token=05982663a539c9fac05608e4b7ca76649356f1cdb3dcc750
Shutdown this notebook server (y/[n])? **y**
[C 18:43:10.271 NotebookApp] Shutdown confirmed
[I 11:48:03.637 NotebookApp] Shutting down 1 kernel
[I 11:48:03.838 NotebookApp] Kernel shutdown: 1b0c8fe9-a2de-45bc-ab5a-bea16e60dc0f

[user@cn0859 ddm]$ **exit**
exit
salloc.exe: Relinquishing job allocation 1151805


```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. ampl.sh). For example:



```

#!/bin/bash
set -e
module load ampl
# command here

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] ampl.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. ampl.swarm). For example:



```

cmd1
cmd2
cmd3
cmd4

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f ampl.swarm [-g #] [-t #] --module ampl
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module ampl Loads the ampl module for each subjob in the swarm 
 | |
 | |
 | |












