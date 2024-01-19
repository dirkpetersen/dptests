

document.querySelector('title').textContent = ' Amber on Biowulf';

Amber on Biowulf



|  |
| --- |
| 
Quick Links
[Batch job on Biowulf](#batch)
[On a single GPU](#gpu)
[On two GPUs](#2gpu)
[Walltime Limits and Chaining](#chain)
[Benchmarks](amber)
[Documentation](#doc)
 |


**[AMBER (Assisted Model Building with Energy
Refinement)](http://ambermd.org)** is a package of molecular simulation programs. AMBER
contains a large number of of modules; note that only the sander modules and
pmemd are parallelized.



The term "Amber" refers to two things. First, it is a set of molecular mechanical force fields for the simulation of biomolecules (these force fields are in the public domain, and are used in a variety of simulation programs). Second, it is a package of molecular simulation programs which includes source code and demos.

Amber is developed by: David Case at Rutgers University, Tom Cheatham at the University of Utah, Ken Merz at Michigan State University, Adrian Roitberg at the University of Florida, Carlos Simmerling at SUNY-Stony Brook, Scott LeGrand at NVIDIA, Darrin York at Rutgers University, Ray Luo at UC Irvine, Junmei Wang at the University of Pittsburgh, Maria Nagan at Stony Brook, Ross Walker at GSK, and many others. Amber was originally developed under the leadership of Peter Kollman.

[Amber website](https://ambermd.org)

[A bug in parmed was reported in Jan 2023.](https://github.com/ParmEd/ParmEd/issues/1280) The patch
will be included in the upcoming AmberTools23. In the meantime, any users who use Amber to build CHARMM forcefields should
use the independent [parmed](parmed.html), which has been patched. For more information, please contact Alex Sodt alexander.sodt@nih.gov


There are several versions of Amber available. The module names are as follows: 


|  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| **Module Name** **Details**| **amber/22.gcc** built with gcc7.4.0, OpenMPI 4.0.4 | |
| **amber/22-gpu** | built with CUDA 11.3, gcc 7.4.0, openmpi/4.0.5. Will run on the A100 GPUs. You may see errors like Note: The following floating-point exceptions are signalling: IEEE\_UNDERFLOW\_FLAG IEEE\_DENORMAL but these can be ignored as per 
<http://archive.ambermd.org/201804/0130.html> |
| **amber/20.intel** built with Intel 2020.0.166 compiler, OpenMPI 4.0.4. This version has better performance than the gcc build (amber/20.gcc), and runs on all of the nodes in norm and multinode, but will not run on the oldest nodes (x2670 and x5660) in the quick queue. 
| **amber/20-gpu** built with gcc 7.4.0, OpenMPI 4.0.4. Intended for jobs that are run on the GPUs. Will run on all GPU types except A100.
| **amber/20.gcc** This is the same build as amber/20-gpu, but does not include the CUDA libraries so will not run on GPUs. It is appropriate for running on norm and multinode if for some reason the default Intel build is not desired. 

| **amber/18** built with the AVX2 instruction set, Intel 2017.4.196 compilers, OpenMPI 2.1.2. This version has better performance, and runs on all of the nodes in norm and multinode, but will not 
run on the oldest nodes (x2670 and x5660) in the quick queue. 
| **amber/18-gpu**  module sets up CUDA library paths as well. Intended for GPU runs.
| **amber/18.generic**  built with Intel 2017.4.196 compilers, OpenMPI 2.1.2. This version is slower because it does not use AVX2, but will run on any node in the 
cluster including the oldest quick queue nodes.
| **amber/16**  built with the AVX2 instruction set, Intel 2017.4.196 compilers, OpenMPI 2.1.2. This version has better performance, and runs on all of the nodes in norm and multinode, but will not 
run on the oldest nodes (x2670 and x5660) in the quick queue. 
| **amber/16-gpu** module sets up CUDA library paths as well. Intended for GPU runs.
| **amber/16.generic**  built with Intel 2017.4.196 compilers, OpenMPI 2.1.2. This version is slower because it does not use AVX2, but will run on any node in the 
cluster including the oldest quick queue nodes.
 | |
 | |
 | |
 | |
 | |
 | |
 | |
 | |
 | |
 | |


Thus, if you wish to run on any node in the quick queue, use the module amber/\*.generic

Setting up Amber files with LEaP
LEaP is a graphical builder of input files for AMBER modules. LEaP can be
used via the Xwindows graphical interface xleap, or the terminal
version tleap. To run xleap,


1. Open an Xwindows session to Biowulf ([More information](/docs/connect.html) about Xwindows on Macs,
Windows, and Unix desktop machines.)
2. Load the module for the version you want, and then type 'xleap'.

```

biowulf% module load amber/16
biowulf% xleap

```

You should see the xleap window appear, in which you can type any LEaP
commands.



See the [AMBER tutorials](http://ambermd.org/tutorials/) for more information.

Batch job on Biowulf

For basic information about setting up an Amber job, see [the Amber manual](http://www.lulu.com/content/2369585) and the [Amber tutorials](http://ambermd.org/tutorials/) . 

The Amber executables can run in parallel on all Biowulf computational nodes. However, 
[benchmark runs](/apps/amber/) indicate that Amber jobs scale best to the CPUs on a 
single node. Therefore we recommend that users run Amber jobs on the regular norm partition nodes or [on the GPU nodes](#gpu). 
To determine the most appropriate number of CPUs to allocate, you should run your own benchmarks. 


**Sample script**


```

#!/bin/bash
# This file is amber.run
#

module load amber/16

cd /data/$USER/amber/myproject
mpirun  $AMBERHOME/bin/pmemd.MPI -O -i mdin -o mdout -inf mdinfo -x mdcrd -r restrt


```

Submit with, for example:


```

sbatch --ntasks=8 --ntasks-per-core=1 --nodes=1 --time=168:00:00 --exclusive amber.run

```

This job would be run on 8 cores of a single node, and will not utilize hyperthreaded cores.
The max walltime is set to 168 hrs, which is a week. See the [section on
walltime limits below](#chain).

On a single GPU

Amber runs extremely fast on a single GPU. Since the GPU performance is significantly better than the CPU performance, it is
worth running most Amber jobs on a single GPU. (see [benchmarks](/apps/amber/index.html)). Larger molecular systems may 
benefit from running on more than 1 GPU, but please run your own benchmarks to make sure (and send them to us!)

Set up your Amber batch script along the following lines:

```

#!/bin/bash

cd /data/$USER/mydir

module load amber/16-gpu

$AMBERHOME/bin/pmemd.cuda -O -i mdin -o mdout -inf mdinfo -x mdcrd -r restrt

```


Submit this job with:

```

sbatch --partition=gpu --gres=gpu:k80:1 jobscript       (1 K80 GPU)
or
sbatch --partition=gpu --gres=gpu:p100:1 jobscript      (1 P100 GPU)
or
sbatch --partition=gpu --gres=gpu:v100:1 jobscript      (1 V100 GPU)

```

where


|  |  |  |  |
| --- | --- | --- | --- |
| --partition=gpu  submit to the GPU partition
|  --gres=gpu:k20x:1  allocate a single k20x GPU for this job
 | |
 | |



The jobload command will show 1 CPU being used. The output from Amber will indicate the GPU usage. The 'nvidia-smi' command can also be used to check whether Amber executables are using the GPU (as described in the section below)

On 2 GPUs

 It is not possible to run a single Amber job on both the K20s on a node, since those 2 GPUs do not have peer-to-peer communication. (see the [Amber GPU page](http://ambermd.org/gpus/recommended_hardware.htm) for an explanation of peer-to-peer communication). 

However, on the K80 nodes, the GPUs do have peer-to-peer communication. It is therefore possible to run on 2 GPUs on a K80 node. However, the performance in most cases is worse on 2 GPUs than on a single GPU. If you plan to run a job on 2 GPUs, please run benchmarks first and verify that the performance is better on 
2 GPUs than on 1. ([Benchmarks](/apps/amber/)). Note that the batch system will set the variable $CUDA\_VISIBLE\_DEVICES to the allocated GPUs.

Sample batch script:

```

#!/bin/bash

module load amber/16-gpu
cd /path/to/your/dir
mpirun -np 2 $AMBERHOME/bin/pmemd.cuda.MPI -O -i in1 -o out1 -inf info1 -x crd1 -r r1 

```

Submit with:

```

sbatch --partition=gpu --gres=gpu:k80:2  --time=12:00:00  jobscript

```




You can check the behaviour of your job with the 'nvidia-smi' utility. 
Determine the GPU node on which your job is running via [jobload](https://hpc.nih.gov/docs/biowulf_tools.html#jobload) or [sjobs](https://hpc.nih.gov/docs/biowulf_tools.html#sjobs). 
Suppose your job is on node cn0626, and is using 2 GPUs:

```

biowulf% **ssh cn0626 nvidia-smi**
Mon Jan 16 17:52:13 2017
+-----------------------------------------------------------------------------+
| NVIDIA-SMI 367.48                 Driver Version: 367.48                    |
|-------------------------------+----------------------+----------------------+
| GPU  Name        Persistence-M| Bus-Id        Disp.A | Volatile Uncorr. ECC |
| Fan  Temp  Perf  Pwr:Usage/Cap|         Memory-Usage | GPU-Util  Compute M. |
|===============================+======================+======================|
|   0  Tesla K80           On   | 0000:83:00.0     Off |                  Off |
| N/A   41C    P0   115W / 149W |    145MiB / 12205MiB |     75%      Default |
+-------------------------------+----------------------+----------------------+
|   1  Tesla K80           On   | 0000:84:00.0     Off |                  Off |
| N/A   27C    P8    33W / 149W |      0MiB / 12205MiB |      0%      Default |
+-------------------------------+----------------------+----------------------+
|   2  Tesla K80           On   | 0000:8A:00.0     Off |                  Off |
| N/A   70C    P0   128W / 149W |    145MiB / 12205MiB |     80%      Default |
+-------------------------------+----------------------+----------------------+
|   3  Tesla K80           On   | 0000:8B:00.0     Off |                  Off |
| N/A   45C    P8    33W / 149W |      0MiB / 12205MiB |      0%      Default |
+-------------------------------+----------------------+----------------------+

+-----------------------------------------------------------------------------+
| Processes:                                                       GPU Memory |
|  GPU       PID  Type  Process name                               Usage      |
|=============================================================================|
|    0     51187    C   ...cal/apps/amber/amber16/bin/pmemd.cuda.MPI    61MiB |
|    0     51188    C   ...cal/apps/amber/amber16/bin/pmemd.cuda.MPI    79MiB |
|    2     51187    C   ...cal/apps/amber/amber16/bin/pmemd.cuda.MPI    79MiB |
|    2     51188    C   ...cal/apps/amber/amber16/bin/pmemd.cuda.MPI    61MiB |
+-----------------------------------------------------------------------------+

```

The GPU numbers reported by 'nvidia-smi' may not match the GPUs you specified with the 'CUDA\_VISIBLE\_DEVICES' variable.



Walltime limits and chaining jobs

Walltime limits are set on most Biowulf partitions. Type 'batchlim' to see the current walltime limits, or see the [systems
status page](/systems/status/). Note that the default walltime on the norm queue is 4 hrs, but you can extend this to 10 days. 
Amber jobs should be designed to run for a week or so, save a checkpoint file, and submit a new job starting from that checkpoint.

An example batch script is below. This script runs a single simulation, saves a copy of the output files, and then resubmits a new job starting 
from Amber's 'restart' file. 


```

#!/bin/bash
# this file is called amber.run

module load  amber/16
module list

echo "Running on $SLURM_NTASKS corse"



# rename the restart file to the coordinate filename
mv restrt inpcrd

#run sander
mpirun -np $SLURM_NTASKS `which sander.MPI` -O -i mdin -c inpcrd -p prmtop -r restrt -x traj -e ene -o mdout

#keep a copy of the output from this run
 mv mdout  md_$run.out
 mv traj  md_$run.trj
 mv ene  md_$run.ene
 cp restrt  md_$run.rst

# if less than 10 runs have been performed, increase the run number and submit the next job
if (( "$run" < "10" ))
   then
     run=`expr $run + 1`
     sbatch --ntasks=8 --ntasks-per-core=1 --time=168:00:00 --exclusive amber.run
fi

```


To submit this job, copy the original input coordinate file to 'restrt' for the first run, and then submit.

```

cp inpcrd restrt               
sbatch --ntasks=8 --ntasks--per-core=1 --time=168:00:00 --exclusive amber.run

```


Benchmarks

Based on the [benchmarks](/apps/amber/index.html), it is highly recomended that you run Amber on a GPU node. 

[Full benchmark details](/apps/amber/)

Documentation

[Amber 20 reference manual](http://ambermd.org/doc12/Amber20.pdf)  

[Amber 18 reference manual](http://ambermd.org/doc12/Amber18.pdf)  

[Amber 16 reference manual](http://ambermd.org/doc12/Amber16.pdf)  

[Amber website](http://ambermd.org)




























































