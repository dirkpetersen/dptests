

document.querySelector('title').textContent = ' NAMD on Biowulf ';

NAMD on Biowulf 



|  |
| --- |
| 
Quick Links
[Infiniband batch job](#parallel)
[On GPUs](#gpu)
[Replica Exchange](#replica)
[NAMD with Plumed](#plumed)
[Walltimes and chaining jobs](#chain)
[Benchmarks](namd/index.html)
[Documentation](#doc)
 |


[NAMD](http://www.ks.uiuc.edu/Research/namd/) is a parallel
molecular dynamics program for UNIX platforms designed
for high-performance simulations in structural biology. It is developed by the
[Theoretical Biophysics Group at the Beckman Center, University of Illinois](https://www.ks.uiuc.edu/).


NAMD was developed to be compatible with existing molecular dynamics
packages, especially the packages X-PLOR and [CHARMM](/apps/charmm/), so it will accept X-PLOR and CHARMM input
files. The output files produced by NAMD are also compatible with X-PLOR and
CHARMM.


NAMD is closely integrated with with [VMD](VMD.html) for
visualization and analysis.



There are several versions of NAMD available on Biowulf. You can check the available versions with

```

module avail NAMD

```

**Important:** Please read the webpage [Making efficient use of Biowulf's Multinode Partition](https://hpc.nih.gov/policies/multinode.html) before running large parallel jobs. 


On Helix

NAMD is a compute-intensive scientific program that cannot be run on Helix. Helix is a dedicated interactive data transfer node and should be used only for file transfer and file management..


Batch job on Biowulf

The following example uses the ApoA1 benchmark example from the NAMD site. It is available on Biowulf in   

/usr/local/apps/NAMD/TESTDATA

**Specifying a homogenous set of nodes**

The 'multinode' partition, to which all jobs that require more than a single node must be submitted, is heterogenous. For efficient multinode 
parallel jobs, you need to ensure that you request nodes of a 
single CPU type. For example, at the time of writing this webpage, the 'freen' command displays:


```

biowulf% freen

Partition    FreeNds    FreeCPUs      Cores  CPUs    Mem   Disk   Features
-------------------------------------------------------------------------------------------------------
...
multinode   65/466     3640/26096        28    56    248g   400g   cpu56,core28,g256,ssd400,x2695,ibfdr
multinode   4/190      128/6080          16    32     60g   800g   cpu32,core16,g64,ssd800,x2650,ibfdr
multinode   312/539    17646/30184       28    56    250g   800g   cpu56,core28,g256,ssd800,x2680,ibfdr
...

```


These lines indicate that there are 3 kinds of nodes in the multinode partition. You should submit your job exclusively to one kind of node by specifying --constraint=x2695, 
--constraint=x2650, or --constraint=x2680 as in the examples below.


Sample batch script for the ibverbs version: 

```

#!/bin/bash

cd /data/$USER/mydir

module load NAMD/2.14-verbs
make-namd-nodelist
charmrun ++nodelist ~/namd.$SLURM_JOBID ++p $SLURM_NTASKS `which namd2` +setcpuaffinity  input.namd

# delete the NAMD-specific node list
rm ~/namd.$SLURM_JOBID

```

**Note:** The NAMD +setcpuaffinity flag should be used for the ibverbs version for performance improvement. This flag should **not** be
used when running the OpenMPI/Intel compiled version, since OpenMPI enforces its own cpu affinity. It should also not be used when you are not allocating all the CPUs on a node, since it assigns the cpu affinity in a round-robin fashion. See [https://www.ks.uiuc.edu/Research/namd/2.9/ug/node87.html](https://www.ks.uiuc.edu/Research/namd/2.9/ug/node87.html
).

Sample batch script for the OpenMPI 2.0/Intel-compiler version compiled on Biowulf   

Note: in our benchmarks, this version was slightly slower than the ibverbs version, so most users will want to use the ibverbs version.

```

#!/bin/bash
#
cd /data/$USER/mydir
module load NAMD/2.14-openmpi

mpirun -np $SLURM_NTASKS `which namd2`    input.namd

```


Submit this job with: 

```

sbatch --partition=multinode --constraint=x2680 --ntasks=# --ntasks-per-core=1 --time=168:00:00 --exclusive jobscript

```


where:

--partition=multinode 
Submit to the IB-FDR partition. Highly parallel programs like NAMD with lots of interprocess communication should be run on the Infiniband (IB) network
--constraint=x2680
Request only x2680 nodes in the multinode partition. A heterogenous set of nodes is likely to lower performance.
--ntasks=# 
Specifies the number of NAMD processes to run. This is passed into the script via the $SLURM\_NTASKS variable.   

The number of tasks should be (number of nodes) \* (number of physical cores per node).   

e.g. for the x2695 nodes, ntasks should be (number of nodes) \* 28.

--ntasks-per-core=1 
Specifies that each NAMD process should run on a physical core. Hyperthreading is ignored. This parameter is highly recommended for parallel jobs. 
--time=168:00:00
Specifies a walltime limit of 168 hrs = 1 week. See the section on [chaining jobs](#chain) below.
--exclusive 
specifies that all nodes allocated to this job should be allocated exclusively. This is recommended for multi-node parallel jobs.


Due to a technical complication, 'jobload' may report incorrect results for a NAMD parallel job. Here is a typical NAMD ibverbs run with the jobload showing as 0:

```

           JOBID            TIME            NODES  CPUS  THREADS   LOAD       MEMORY
                     Elapsed / Wall               Alloc   Active           Used /     Alloc
        32072214    00:03:29 /    08:00:00 cn1517    56        0     0%     0.0 /   56.0 GB
                    00:03:29 /    08:00:00 cn1518     0        0     0%     0.0 /    0.0 GB
                 Nodes:    2    CPUs:  112  Load Avg:   0%

```

However, 'ssh cn1517 ps -C namd2' will show that there are 28 namd2 processes on each node. The NAMD output file will also report details such as:

```

Charm++> Running on 2 unique compute nodes (56-way SMP).

```


On GPUs

The latest GPU-resident single-node-per-replicate computation mode NAMD 3.0 build is available via the module 
NAMD/3.0\*-CUDA. Please read
the [NAMD webpage about NAMD 3.0](https://www.ks.uiuc.edu/Research/namd/3.0/announce.html)
Note that to get the enhanced performance, you
must run this version on a single node with NVLink, which means the v100x or a100 GPU nodes on Biowulf. 

**Sample script, single-node GPU job, NAMD 3.0beta3**

```

#!/bin/bash

cd /data/$USER/mydir
module load NAMD/3.0/beta3-CUDA
namd3 +p${SLURM_CPUS_PER_TASK} +setcpuaffinity stmv.namd

```

Submit with, for example: 

```

sbatch -p gpu --gres=gpu:v100x:1 --cpus-per-task=8 my_batch_script	# 1 v100x GPU
sbatch -p gpu --gres=gpu:a100:2 --cpus-per-task=16 --nodes=1 my_batch_script	# 2 a100 GPU
sbatch -p gpu --gres=gpu:v100x:4 --cpus-per-task=32 --nodes=1 my_batch_script	# all 4 v100x GPU on a single node

```


See the [benchmarks page](namd) for some samples of performance.

**Single-node GPU job, NAMD 2.14**   

To run a single-node GPU job, that will run on a single K80, P100, V100 or V100x node, create a batch script along the 
following lines. Note thatn NAMD 2.14 can run (and it is recommended to) multiple ntasks per GPU.

```

#!/bin/bash

cd /data/$USER/mydir

module load NAMD/2.14-verbs-CUDA
charmrun ++local `which namd3` +p $SLURM_NTASKS +setcpuaffinity stmv.namd

```

The environment variable $CUDA\_VISIBLE\_DEVICES will be set by Slurm to the GPU devices that are allocated to the job.

To submit to 2 GPU devices and half the CPUs on a K80 node:

```

sbatch --partition=gpu --gres=gpu:k80:2 --ntasks=14 --ntasks-per-core=1 jobscript

```


To submit to all 4 GPU devices and all the CPUs on a V100 node:

```

sbatch --partition=gpu --gres=gpu:v100:4 --ntasks=28 --ntasks-per-core=1 --exclusive jobscript

```


As per the [NAMD GPU documentation](http://www.ks.uiuc.edu/Research/namd/2.12/ug/node90.html), multiple NAMD threads can utilize the same set of GPUs, and the tasks 
are equally distributed among the allocated GPUs on a node. 

**Multi-node GPU job**


While it is possible to run a multinode GPU NAMD job, please be sure that your NAMD job scales to more than 1 GPU node before submitting multinode GPU jobs. ([See 
our benchmarks for details](namd/)).
To submit a multinode job, you could use a script like the following:

```

#!/bin/bash

cd /data/$USER/mydir

module load  NAMD/2.14-verbs-CUDA

# on a K80 node
make-namd-nodelist
charmrun ++nodelist ~/namd.$SLURM_JOBID ++p $SLURM_NTASKS `which namd2` ++ppn 28 input.namd

```

To submit to 2 K80 nodes:

```

sbatch ---partition=gpu --gres=gpu:k80:4 --ntasks=56 --ntasks-per-core=1 --nodes=2 --exclusive jobscript  

```

Note that the number of ntasks is set to the number of cores on 2 nodes, i.e. 56. 

**Monitoring GPU jobs**

To monitor your GPU jobs, use 'jobload' to see the CPU utilization (should be ~ 50%), and 'ssh nodename nvidia-smi' to see the GPU utilization. 
In the example below, a NAMD job is submitted to 8 GPUs (2 nodes) and 56 cores (all cores on the 2 nodes). 

```

[biowulf]$  **sbatch --partition=gpu --gres=gpu:k80:4 --ntasks=56 --ntasks-per-core=1 --nodes=2 run.gpu**
129566

```

Jobload shows that the job is utilizing all cores:

```

[biowulf]$  **jobload -u user**
     JOBID      RUNTIME     NODES   CPUS    AVG CPU%            MEMORY
                                                              Used/Alloc
    129566     00:00:26    cn0603     56       50.00    836.9 MB/62.5 GB
               00:00:26    cn0604     56       50.06    644.6 MB/62.5 GB

```

The 'nvidia-smi' command shows that there are 4 NAMD processes running on the 4 GPUs of the node. The 'GPU-Util' value will bounce around, so is not very meaningful.

```

[biowulf]$  **ssh cn3084 nvidia-smi**
Sun Feb 26 15:19:07 2017
+-----------------------------------------------------------------------------+
| NVIDIA-SMI 367.48                 Driver Version: 367.48                    |
|-------------------------------+----------------------+----------------------+
| GPU  Name        Persistence-M| Bus-Id        Disp.A | Volatile Uncorr. ECC |
| Fan  Temp  Perf  Pwr:Usage/Cap|         Memory-Usage | GPU-Util  Compute M. |
|===============================+======================+======================|
|   0  Tesla K80           On   | 0000:83:00.0     Off |                  Off |
| N/A   48C    P0    58W / 149W |     91MiB / 12205MiB |     15%      Default |
+-------------------------------+----------------------+----------------------+
|   1  Tesla K80           On   | 0000:84:00.0     Off |                  Off |
| N/A   35C    P0    75W / 149W |     91MiB / 12205MiB |     15%      Default |
+-------------------------------+----------------------+----------------------+
|   2  Tesla K80           On   | 0000:8A:00.0     Off |                  Off |
| N/A   51C    P0    62W / 149W |     90MiB / 12205MiB |     12%      Default |
+-------------------------------+----------------------+----------------------+
|   3  Tesla K80           On   | 0000:8B:00.0     Off |                  Off |
| N/A   39C    P0    76W / 149W |     91MiB / 12205MiB |     14%      Default |
+-------------------------------+----------------------+----------------------+

+-----------------------------------------------------------------------------+
| Processes:                                                       GPU Memory |
|  GPU       PID  Type  Process name                               Usage      |
|=============================================================================|
|    0     40816    C   ..._2.12_Linux-x86_64-ibverbs-smp-CUDA/namd2    87MiB |
|    1     40816    C   ..._2.12_Linux-x86_64-ibverbs-smp-CUDA/namd2    87MiB |
|    2     40816    C   ..._2.12_Linux-x86_64-ibverbs-smp-CUDA/namd2    86MiB |
|    3     40816    C   ..._2.12_Linux-x86_64-ibverbs-smp-CUDA/namd2    87MiB |
+-----------------------------------------------------------------------------+

```

Replica Exchange

**Sample Replica Exchange job script using the OpenMPI version on Infiniband (multinode) partition**

```

#!/bin/bash

cd /data/$USER/mydir
module load NAMD/2.13b2-openmpi

mkdir output
(cd output; mkdir 0 1 2 3 4 5 6 7)
mpirun namd2 +replicas 8 job0.conf +stdout output/%d/job0.%d.log

```

The number of MPI ranks must be a multiple of the number of replicas. Thus, for the 8 replicas above, you could submit with:

```

sbatch --partition=multinode --ntasks=24 --ntasks-per-core=1 --nodes=1 --exclusive   jobscript

```

using 24 of the 28 physical cores on a single node. 

**Replica Exchange job using the verbs-CUDA version on multiple GPUs**   


```

#!/bin/bash

cd $SLURM_SUBMIT_DIR

# this run uses the fold_alanin example provided by NAMD
# download the example files and set up the output directories
tar xvf /usr/local/apps/NAMD/2.13b2/replica_example.tar.gz
cd replica/example
mkdir output
(cd output; mkdir 0 1 2 3 4 5 6 7)

module load NAMD/NAMD/2.13b2-verbs-CUDA
make-namd-nodelist

charmrun ++nodelist ~/namd.$SLURM_JOBID +p8 \
    `which namd2` +replicas 8 job0.conf +stdout output/%d/job0.%d.log

```

To run the above on 8 k80 GPUs (2 nodes), you would submit with: 

```

sbatch --partition=gpu --gres=gpu:k80:4 --nodes=2 --ntasks=32 --exclusive  jobscript

```

Note that jobload will report incorrect usage for this job. It will look like: 

```

           JOBID            TIME            NODES  CPUS  THREADS   LOAD       MEMORY
                     Elapsed / Wall               Alloc   Active           Used /     Alloc
        48970813    00:16:44 /    02:00:00 cn4200    56        0     0%     0.0 /  112.0 GB
                    00:16:44 /    02:00:00 cn4201     0        0     0%     0.0 /    0.0 GB

```

However, the appropriate processes and GPU usage can be checked with commands such as the following. For the example above (8 replicas, 8 GPUs), you should see 
4 namd2 processes on the node CPUs, and 4 namd2 processes on the GPUs of each node.

```

biowulf% **ssh cn4200 ps auxw | grep namd**
user   51222  0.0  0.0  17168  1496 ?        S    11:18   0:00 charmrun ++nodelist /home/user/namd.48970813 +p8 /usr/local/apps/NAMD/2.13b2/NAMD_2.13_Linux-x86_64-verbs-smp-CUDA/namd2 +replicas 8 job0.conf +stdout output/%d/job0.%d.log
user   51297  100  0.1 323977568 299836 ?    Rl   11:18   4:29 /usr/local/apps/NAMD/2.13b2/NAMD_2.13_Linux-x86_64-verbs-smp-CUDA/namd2 +replicas 8 job0.conf +stdout output/%d/job0.%d.log
user   51299  100  0.1 323977572 297048 ?    Rl   11:18   4:29 /usr/local/apps/NAMD/2.13b2/NAMD_2.13_Linux-x86_64-verbs-smp-CUDA/namd2 +replicas 8 job0.conf +stdout output/%d/job0.%d.log
user   51300  100  0.1 323977572 297344 ?    Rl   11:18   4:29 /usr/local/apps/NAMD/2.13b2/NAMD_2.13_Linux-x86_64-verbs-smp-CUDA/namd2 +replicas 8 job0.conf +stdout output/%d/job0.%d.log
user   51301  100  0.1 323977568 296996 ?    Rl   11:18   4:29 /usr/local/apps/NAMD/2.13b2/NAMD_2.13_Linux-x86_64-verbs-smp-CUDA/namd2 +replicas 8 job0.conf +stdout output/%d/job0.%d.log

biowulf% **ssh cn4200 nvidia-smi**
Wed Feb 19 11:23:33 2020
+-----------------------------------------------------------------------------+
| NVIDIA-SMI 440.33.01    Driver Version: 440.33.01    CUDA Version: 10.2     |
|-------------------------------+----------------------+----------------------+
| GPU  Name        Persistence-M| Bus-Id        Disp.A | Volatile Uncorr. ECC |
| Fan  Temp  Perf  Pwr:Usage/Cap|         Memory-Usage | GPU-Util  Compute M. |
|===============================+======================+======================|
|   0  Tesla K80           On   | 00000000:83:00.0 Off |                  Off |
| N/A   64C    P0    75W / 149W |    414MiB / 12206MiB |     91%      Default |
+-------------------------------+----------------------+----------------------+
|   1  Tesla K80           On   | 00000000:84:00.0 Off |                  Off |
| N/A   30C    P8    34W / 149W |     11MiB / 12206MiB |      0%      Default |
+-------------------------------+----------------------+----------------------+
|   2  Tesla K80           On   | 00000000:8A:00.0 Off |                  Off |
| N/A   42C    P8    27W / 149W |     11MiB / 12206MiB |      0%      Default |
+-------------------------------+----------------------+----------------------+
|   3  Tesla K80           On   | 00000000:8B:00.0 Off |                  Off |
| N/A   37C    P8    34W / 149W |     11MiB / 12206MiB |      0%      Default |
+-------------------------------+----------------------+----------------------+

+-----------------------------------------------------------------------------+
| Processes:                                                       GPU Memory |
|  GPU       PID   Type   Process name                             Usage      |
|=============================================================================|
|    0     51297      C   ..._2.13_Linux-x86_64-verbs-smp-CUDA/namd2   100MiB |
|    1     51299      C   ..._2.13_Linux-x86_64-verbs-smp-CUDA/namd2   100MiB |
|    2     51300      C   ..._2.13_Linux-x86_64-verbs-smp-CUDA/namd2   100MiB |
|    3     51301      C   ..._2.13_Linux-x86_64-verbs-smp-CUDA/namd2   100MiB |
+-----------------------------------------------------------------------------+

```


NAMD with plumed
NAMD version 2.9 has been compiled with support for [Plumed 2.4.2](http://www.plumed.org/), a library for performing
free energy calculation as part of molecular simulations. To use Plumed, you must load the NAMD/2.9-plumed module. This
version is compiled with OpenMPI to allow parallelization of Plumed, and therefore, the sample batch job for the OpenMPI version
given [above](#parallel) should be used as the basis of scripts.


An example submission script using Plumed would look like:




```

#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=16
#SBATCH --ntasks=64
#SBATCH --partition=multinode
#SBATCH --constraint=x2650
#SBATCH --mem=48g

module load NAMD/2.9-plumed

mpirun -np $SLURM_NTASKS --mca btl self,sm,openib --mca btl_openib_if_exclude "mlx4_0:2" `which namd2` input 

```

Note that the necessary compiler, MPI, and FFTW modules are loaded automatically by the NAMD module.



Walltimes and chaining jobs

There are walltime limits on most Biowulf partitions. Use 'batchlim' to see the current walltime limits. 

 An example namd config file for running a second simulation starting from the last timestep and the restart files of a previous simulation is available at <http://www.ks.uiuc.edu/~timisgro/sample.conf>.

If restarting a NAMD REMD job, be sure to comment out the 'bincoordinates' and 'extendedsystem' parameters in your NAMD configuration file, if applicable

After an initial run has produced a set of restart files, you would submit future runs using a batch script along these lines:


```

#!/bin/bash

module load NAMD/2.10

# Create host file (required)
make-namd-nodelist
mpirun -n $SLURM_NTASKS  `which namd2` myjob.restart.namd > out.log
rm -f ~/namd.$SLURM_JOBID

# this script resubmits itself to the batch queue. 
# The NAMD config file is set up to start the simulation from the last timestep 
#   in the previous simulation
sbatch --partition=multinode --constraint=x2650 --ntasks=$SLURM_NTASKS --ntasks-per-core=1 --time=168:00:00 --exclusive this_job_script

```

Submit this script, as usual, with a command like:

```

sbatch --partition=multinode --constraint=x2650 --ntasks=64 --ntasks-per-node=1 --time=168:00:00 --exclusive this_job_script

```

The NAMD 2.10 replica.namd file is at /usr/local/apps/NAMD/NAMD\_2.10\_Linux-x86\_64-ibverbs/lib/replica/replica.namd.

#### Chaining Replica Exchange jobs


(Thanks to Christopher Siwy (CC) for this information).

* Create separate job/config/namd files for running the restart job
* Modify them appropriately (e.g., reference the correct restart files, do not remove the existing replica output folders)
* Update the job0.conf file to source the restart tcl file that NAMD generates in the folder up from each replica's folder



The most important points here are:  

- Ensure you do not delete your replica folders when you run the restart (as this is usually done when you start a new REMD simulation)  

- In your job0.conf (or whatever you name it) file, include the following two lines after the line referencing the NAMD configuration  


```

source [format $output_root.job0.**restart20.tcl** ""]
set num_runs **10000**
 
```

The items in bold will likely be subjective. And the number of runs, is the TOTAL number of runs for the simualtion, not the number of 
runs to run from that point forward. So in the example above, the restart will begin at the 20th run and continue till it reaches the 10,000th run.

[This thread in the NAMD mailing list](http://www.ks.uiuc.edu/Research/namd/mailing_list/namd-l.2013-2014/1719.html) may help in debugging problems.


Benchmarks

[On the Benchmarks page](namd/)

Documentation


[Theoretical and Computational Biophysics group](http://www.ks.uiuc.edu) at UIUC, the NAMD/VMD developers.


























































































