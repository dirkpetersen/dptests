

document.querySelector('title').textContent = ' Gromacs on Biowulf';

Gromacs on Biowulf



|  |
| --- |
| 
Quick Links
[Batch job on Biowulf](#parallel)
[Replica Exchange](#replica)
[on GPUs](#gpu)
[Walltimes & Chaining jobs](#chaining)
[Tips for Best Performance](#tips)
[Benchmarks](/apps/gromacs/)
[Documentation](#doc)
 |


GROMACS ([www.gromacs.org](http://www.gromacs.org)) is a
versatile package to perform molecular dynamics, i.e. simulate the Newtonian
equations of motion for systems with hundreds to millions of particles. It is
primarily designed for biochemical molecules like proteins and lipids that have
a lot of complicated bonded interactions, but since GROMACS is extremely fast
at calculating the nonbonded interactions (that usually dominate simulations)
many groups are also using it for research on non-biological systems, e.g.
polymers.



Important Notes
* Module Name: gromacs (see [the modules page](/apps/modules.html) for more information)
* MPI Parallelization. Please read the webpage [Making efficient use of Biowulf's Multinode Partition](https://hpc.nih.gov/policies/multinode.html) before running large parallel jobs.




Batch job on Biowulf

Gromacs can multi-thread as well as use MPI. For small jobs, e.g. 8 cpus on a single node, multi-threading works almost as well as MPI. For larger jobs, it is best to use MPI. 
(see the [Benchmarks page](gromacs/index.html)) for details. 

**Specifying a homogenous set of nodes**

The 'multinode' partition, to which all jobs that require more than a single node must be submitted, is heterogenous. For efficient parallel jobs, you need to ensure that you request nodes of a 
single CPU type. For example, at the time of writing this webpage, the 'freen' command displays:


```

biowulf% freen
...
multinode   65/466       3640/26096        28    56    248g   400g   cpu56,core28,g256,ssd400,x2695,ibfdr
multinode   4/190        128/6080          16    32     60g   800g   cpu32,core16,g64,ssd800,x2650,ibfdr
multinode   312/539      17646/30184       28    56    250g   800g   cpu56,core28,g256,ssd800,x2680,ibfdr
...

```

These lines indicate that there are 3 kinds of nodes in the multinode partition. You should submit your job exclusively to one kind of node by specifying --constraint=x2695 or --constraint=x2680, as in the examples below. (Note that Gromacs 2018.3 and 2020.2 will not currently run on the x2650 nodes)


**Sample MPI batch script for Gromacs 2022.4**

```

#!/bin/bash

module load gromacs/2022.4

mpirun gmx_mpi mdrun -ntomp 1 -s topol.tpr

```

Gromacs will use GPUs if available, but will run on CPUs if not. Therefore, the same batch script will work for both GPUs and CPUs. Sample submission commands: 


```

sbatch --ntasks=# --ntasks-per-core=1 --nodes=1 run.2022.4  # CPUs

sbatch -p gpu --gres=gpu:p100:1 --ntasks=1 --ntasks-per-core=1 run.2022.4  # for 1 p100 GPU

```

where 'p100' can be replaced by whichever GPU type is desired. Use 'freen' to see GPUs available.


**Sample MPI batch script for Gromacs 2021.3 and before:**

```

#!/bin/bash

module load gromacs/2018

mpirun -np $SLURM_NTASKS `which mdrun_mpi` -ntomp 1 -s ion_channel.tpr -maxh 0.50 \
      -resethway -noconfout -nsteps 1000

```

Submit this job with:

```

sbatch --partition multinode --constraint=x2695 --job-name=gmx  --ntasks=# --ntasks-per-core=1  --time=168:00:00 --exclusive jobscript

```

where


|  |  |  |  |  |  |  |  |  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| --partition multinode  Submit to the multinode partition where all nodes are Infiniband-connected
| --constraint=x2695  All nodes should be x2695's. 
| --ntasks #  the number of MPI processes you wish to run. 
| --ntasks-per-core=1  ensures that Gromacs will only run 1 MPI process per physical core (i.e will not use both hyperthreaded CPUs). This is recommended for parallel jobs.| -ntomp1   uses only one OMP thread per MPI thread. This means that Gromacs will run using only MPI, which provides the best performance. 
| --time=168:00:00  max walltime=168 hrs (1 week). See the section on [chaining jobs](#chaining) below.
| --exclusive  Allocate the nodes exclusively to this job (recommended for parallel jobs)
 | |
 | |
 | |
 | |
 | |
 | |
 | |



Replica Exchange

Sample batch script for Gromacs 4.6.5 (thanks to Mingzhen Zhang)


```

#!/bin/bash

module load gromacs/2018

cd $SLURM_SUBMIT_DIR

mpirun -np $SLURM_NTASKS `which mdrun_mpi` -ntomp 1 -s cmd_.tpr -maxh 0.50 -resethway -noconfout -cpi state.cpt -noappend -multi 48 -replex 1000

```


Submit with: 

```

sbatch --partition multinode --constraint=x2695 --job-name=MyJob --ntasks=64 --ntasks-per-core=1 --exclusive   myjobscript

```


On GPUs

GPU support is built in to Gromacs 5.\*. 
Sample batch script with the adh\_cubic job from the [Gromacs GPU documentation](https://manual.gromacs.org/current/user-guide/mdrun-performance.html#running-mdrun-with-gpus). The files for this job are 
available in /usr/local/apps/gromacs/ADH\_bench\_systems.tar.gz (untar the file and look in the adh\_cubic directory).

Note that the only following versions on Biowulf are built with CUDA/GPU support:
* gromacs/2018.3 
* gromacs/2021.3



Other Gromacs versions will be visible when you run 'module avail gromacs', but those versions were built for a single user or lab and were not built with CUDA/GPU support.


Sample batch script:

```

#!/bin/bash

module load gromacs/2018.3

mkdir /data/$USER/gromacs/
cd /data/$USER/gromacs
tar xvzf /usr/local/apps/gromacs/adh_cubic.tar.gz
cd adh_cubic

mpirun -np $SLURM_NTASKS --mca btl_openib_if_exclude "mlx4_0:1" `which mdrun_mpi` \
       -ntomp $SLURM_CPUS_PER_TASK -s topol.tpr

```



Submit this job with, for example:

```

sbatch --partition=gpu --gres=gpu:k80:2 --ntasks=2  --ntasks-per-core=1 --cpus-per-task=1  --time=HH:MM:SS jobscript

```


The above command will submit the job to 2 GPUs. 



|  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| -ntomp 1  sets 1 CPU thread per MPI (GPU) process. In our tests a 1:1 CPU:GPU ratio gave the best performance (see [benchmarks](gromacs/index.html#gpu)).
| --mca btl\_openib\_if\_exclude "mlx4\_0:1" prevents a warning about OpenFabrics from appearing in your output. You can also leave it out
and live with the warning :-).
| --partition=gpu  submit to the GPU partition
| --gres=gpu:k20x:2  Resource=gpu, Resource type=k20x, Count=2 (2 GPUs). Note that 'count' is required even if you are submitting to a single GPU.
| --ntasks=2 Number of MPI tasks to spawn.
| --ntasks-per-core=1  Run only 1 MPI task per physical core. This is best for most MD jobs
| --time=HH:MM:SS  walltime to be allocated for the job -- HH hours, MM minutes, SS seconds
| --cpus-per-task=1 >Number of threads per MPI task. You should run your own benchmarks to determine the best values for ntasks and cpus-per-task. 
 | |
 | |
 | |
 | |
 | |
 | |
 | |
 | |



Chaining jobs

The max walltime on the multinode partition is 10 days. (type 'batchlim' to see the CPU and walltime limits on all partitions).
Thus, jobs should be designed to run for a week or so, save a checkpoint file, and submit a new job starting from that checkpoint. 

 A reasonable strategy would be to set up a job to run for a week or less by setting the number of steps appropriately, and then, at the end of the job, have it resubmit itself to continue the simulation. Below is a sample batch script:


```

#!/bin/bash
# this script is called Run.ib

module load gromacs/2018.3

cd /path/to/my/dir

mpirun -np $SLURM_NTASKS `which mdrun_mpi` -ntomp 1 -s ion_channel.tpr -maxh 0.50 -resethway -noconfout -nsteps 1000

# use tpbconv to create a new topol.tpr file with an increased number of steps
tpbconv -s topol.tpr -extend 500 -o topol2.tpr

#move the newly created topol.tpr into place
mv topol.tpr topol.tpr.prev; mv topol2.tpr topol.tpr

#resubmit this script
sbatch --partition multinode --constraint=x2680 --job-name=gmx  --ntasks=# --ntasks-per-core=1  --time=168:00:00 --exclusive  Run.ib

```

More information at [Extending Simulations](http://www.gromacs.org/Documentation/How-tos/Extending_Simulations) on the Gromacs site.

 If a Gromacs job is terminated unexpectedly (for example, the walltime limit was hit before the mdrun completed), it is simple to restart. The state.cpt file contains all the information necessary to continue the simulation. Use the '-cpi' and '-append' options to mdrun, which will append to existing energy, trajectory and log files. For example:

```

mpirun -n $np `which mdrun_mpi` -s topol.tpr -cpi state.cpt -append

```


More information at [Doing Restarts](http://www.gromacs.org/Documentation/How-tos/Doing_Restarts) on the Gromacs website.

Tips for Best Performance

* Make sure you request homogenous resources. The 'freen' command will show several kinds of nodes in the
multinode partition. Pick one (depending on CPU speed or availability), and submit to only that type of node by using the 
'constraint' flag. e.g.

```

#SBATCH --constraint=x2695

```

Running on a mix of node types will effectively mean running on the slowest type of node. 

* The bottleneck for MD MPI-based programs is the inter-node communication, so you should submit to as few nodes
as possible, utilizing all the cores on the node. For example, suppose you want to run 100 ntasks. 'freen' shows nodes with 
16 or 28 cores. If you run on 28-core nodes, you can utilize 4 full nodes with:

```

#SBATCH --ntasks=112
#SBATCH --ntasks-per-core=1
#SBATCH --nodes=4
#SBATCH --exclusive
#SBATCH --constraint=nodetype

```

* Check the job utilization with 'jobload' -- ideally you should see all the allocated cores fully utilized. Since only one 
process is running per core, you should see 50% utilization on each of the allocated nodes. 

```

**Example of a well set up job**

biowulf% jobload -j 11111111
           JOBID            TIME            NODES  CPUS  THREADS   LOAD       MEMORY
                     Elapsed / Wall               Alloc   Active           Used /     Alloc
      11111111    0-04:05:52 /  2-02:00:00 cn1135    56       28   50%     0.8 /    4.0 GB
                  0-04:05:52 /  2-02:00:00 cn2170    56       28   50%     0.8 /    4.0 GB
                  0-04:05:52 /  2-02:00:00 cn2171    56       28   50%     0.8 /    4.0 GB
                  0-04:05:52 /  2-02:00:00 cn2172    56       28   50%     0.8 /    4.0 GB
  
```


```

**Example of a poorly set up job**

biowulf% jobload -j 11111111
           JOBID            TIME            NODES  CPUS  THREADS   LOAD       MEMORY
                     Elapsed / Wall               Alloc   Active           Used /     Alloc
      11111111    0-04:05:52 /  2-02:00:00 cn1135     4        4   100%     0.3 /    4.0 GB
                  0-04:05:52 /  2-02:00:00 cn2170    16       13    81%     0.8 /    4.0 GB
                  0-04:05:52 /  2-02:00:00 cn2171    16       13    81%     0.8 /    4.0 GB
                  0-04:05:52 /  2-02:00:00 cn2172    48       13    27%     0.8 /    4.0 GB
                  0-04:05:52 /  2-02:00:00 cn3024    44       13    30%     0.8 /    4.0 GB
                  0-04:05:52 /  2-02:00:00 cn3025    44       13    30%     0.8 /    4.0 GB
 
```


Benchmarks

see the [benchmarks page](gromacs/index.html)
Documentation

[Gromacs website](http://www.gromacs.org)  

[mdrun documentation for v5.0.4](http://manual.gromacs.org/programs/gmx-mdrun.html)
































































