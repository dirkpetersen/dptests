

document.querySelector('title').textContent = 'paraview on Biowulf';
paraview on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Client/Server mode](#clientserver) 
[Batch mode](#sbatch)
 |



ParaView is an open-source, multi-platform data analysis and visualization application. It
can be run in different modes:



* In Destop mode, Paraview is run locally on your machine.
* In client/server mode, a ParaView server runs on the biowulf cluster and
 an interactive ParaView client on your computer connects to the server
 for analysis via an ssh tunnel.
* In batch mode, ParaView is used as a framework for analysis and visualisation 
 in biowulf batch jobs through its Python interface.



 ParaView can be run on k80 GPU nodes (egl offscreen rendering) or, without GPU acceleration, on multinode
 (mesa offscreen rendering).

#### References:


* J. Ahrens, B. Geveci, C. Law. *ParaView: An End-User Tool for Large Data Visualization*.
 Visualization Handbook, 2005.


Documentation
* ParaView Main Site: [paraview.org](https://www.paraview.org)
* ParaView Tutorial: [paraview.org](https://www.paraview.org/Wiki/images/1/13/ParaViewTutorial52.pdf)


Important Notes
* Module Name: paraview (see [the modules page](/apps/modules.html) for more information)
* ParaView is an MPI program. Current versions are precompiled and include their own custom MPICH
* For client/server mode please download ParaView from [paraview.org/download](https://www.paraview.org/download/). The version has to match exactly the server version running on biowulf.



Client/Server mode
#### k80 GPU


In this example we will run a ParaView server as a batch job allocating one
node with 2 K80 GPUs each used for hardware accelerated rendering. We then
connect a ParaView client running on the desktop to the server for interactive
analysis. Paraview can either use half a node (1GPU + CPUS), or multiples
of whole nodes (i.e. 1, 2, ... nodes allocated exclusively using both GPUs).
That is due to the way GPUs are utilized by ParaView. The batch example
uses half a node to illustrate how to use $CUDA\_VISIBLE\_DEVICES. The 1 node
client/server example can easily be extended to more than one node. 


To start a ParaView server, create a batch script similar to the following:

```

#! /bin/bash
# this file is paraview_server.sh
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --ntasks-per-core=1
#SBATCH --partition=gpu
#SBATCH --gres=gpu:k80:2
#SBATCH --exclusive
#SBATCH --mem=120g
set -e


module load paraview/5.10.1 || exit 1
pvopts="--disable-xdisplay-test --force-offscreen-rendering"

ntasks_per_gpu=$(( SLURM_NTASKS_PER_NODE / 2 ))

mpiexec -map-by node -iface enet0 \
    -np ${ntasks_per_gpu} pvserver ${pvopts} --displays=0 : \
    -np ${ntasks_per_gpu} pvserver ${pvopts} --displays=1


```

Note that the mpiexec is set up such that half the tasks are assigned to each of the
two GPUs on each node. This will only work properly if only whole nodes are allocated.


In this example we will allocate just 1 node with 2GPUs. 8 tasks each share
a single GPU for rendering. In some circumstances it may be better to use
2 tasks per core, for example. Since we are using sbatch directives we can submit with



```

[user@biowulf]$ **sbatch paraview\_server.sh**
50077458

```

#### Multinode CPU rendering


In some cases the server may need more memory than available on the k80 nodes
or exceed the k80 GPU memory. In such cases the mesa based CPU server can be used.



```

#! /bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=28
#SBATCH --partition=multinode
#SBATCH --exclusive
#SBATCH --mem=240g
#SBATCH --constraint=x2695

set -e

module load paraview/5.10.1-mesa

pvopts="--disable-xdisplay-test --force-offscreen-rendering"
mpiexec -iface ib0 pvserver ${pvopts}

```

#### Client connections


The server will run until the job hits its time limit, is canceled, or is
stopped when the client terminates its connection to the server.


Wait until the server job starts and check the output file for the
port and hostname of the server



```

[user@biowulf]$ **squeue -j 50077458**
   JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
50077458       gpu paraview     user PD       0:00      1 (None)

[user@biowulf]$ **squeue -j 50077458**
   JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
50077458       gpu paraview     user  R       1:36      1 cn0605

[user@biowulf]$ **cat slurm-50077458.out**
[+] Loading paraview 5.10.1
Waiting for client...
Connection URL: cs://cn0605:11111
Accepting connection(s): cn0605:11111

```

Now we need a tunnel from your local machine to the compute node
indicated by the pvserver process. Assuming your machine is on the NIH
campus or on VPN this command will set up the tunnel on a Mac or
a Linux machine to the node via biowulf:



```

[myComputer]$ **ssh -L 11111:localhost:11111 user@biowulf.nih.gov \
 -t ssh -L 11111:localhost:11111 user@cn0605**

```

Replacing 'user' with your username and 'cn0605' with the node
shown in the slurm output file.



Using putty on a Windows machine involves a 2-step process:


* Establish a connection to biowulf with a tunnel from the putty GUI
* Extend the tunnel from biowulf to the compute node from the biowulf command line.


Once a tunnel has been established, start your local ParaView client and click the
connect button. The brief screencast below shows how to connect and open an example
data set.


![](/images/video-icon.png) [Paraview client connection to server](https://www.youtube.com/watch?v=NYfiUcqZktU) (4 min)

Notes:


* Disconnecting from the server terminates the server job.
* The screencast was recored with the Linux version of ParaView 5.4.1.
 Interfaces on other Systems may differ slightly.
* Note that the video uses some properties used in the video only become
visible when selecting the gear icon (advanced properties).
* If port 11111 is used already you can select a different port with the `--server-port`
 option to pvserver
* For server jobs with more than 1 node please make sure to check the slurm
 log file to identify on which node the process listening for client connections
 is running



Batch mode
We will use a very simple python script that renders a sphere and colors
it by the rank of the MPI process that created it:



```

# this file is sphere.py
from paraview.simple import *

sphere = Sphere()
sphere.ThetaResolution = 32
rep = Show()
ColorBy(rep, ('POINTS', 'vtkProcessId'))
Render()
rep.RescaleTransferFunctionToDataRange(True)
Render()
WriteImage('sphere.png')

```

Then create a batch input file to run the python script on a single K80 GPU:



```

#!/bin/bash
# this file is sphere.sh
module load paraview || exit 1

mpirun --mca btl self,vader -np $SLURM_NTASKS \
    pvbatch shpere.py

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```

[user@biowulf]$ **sbatch --ntasks=4 --partition=gpu --mem=10g \
 --gres=gpu:k80:1 --nodes=1 sphere.sh**

```









