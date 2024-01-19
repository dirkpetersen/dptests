

document.querySelector('title').textContent = 'Psi4 on Biowulf';
Psi4 on Biowulf



|  |
| --- |
| 
Quick Links
[On Helix](#helix)
[Multithreaded job on Biowulf](#thread)
[Documentation](#doc)
 |



Psi4 is an ab-initio electronic structure code that supports various methods
for calculating energies and gradients of molecular systems. It is designed
to be highly flexible and extensible, with performance-critical code implemented
in C++ that can be accessed via either Python language bindings or a domain-specific 
input script language. Specific features include:

* Support for various ab initio and Density Functional Theory (DFT) methods
* Numerous coupled cluster methods implemented
* Broad support for Symmetry Adapted Perturbation (SAPT) methods
* Gradients (first derivatives) and Hessian (second derivatives) implemented - allows for efficient geometry optimization of systems
* Allows exporting frequencies in Molden format for visualization of normal modes
* Shared-memory parallelization to run efficiently on multi-core machines



Psi4 is distributed under the GNU Lesser General Public License (LGPL)
version 3. 

Example input files may be found in /usr/local/apps/psi4/1.1/share/psi4/samples.

On Helix

Psi4 is a cpu-intensive program that is not allowed to be run on Helix.

Multithreaded Batch job on Biowulf


Sample batch script: (this file is called run\_psi4)

```

#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --threads-per-core=1
#SBATCH --gres=lscratch:50
#SBATCH --partition=norm

cd $SLURM_SUBMIT_DIR

module load psi4

psi4 -n $SLURM_CPUS_PER_TASK saptest.in

```

Submit this job with 

```

biowulf% sbatch run_psi4

```


where:


|  |  |  |  |  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| -n $SLURM\_CPUS\_PER\_TASK tells Psi4 to run $SLURM\_CPUS\_PER\_TASK threads. This is only defined if --cpus-per-task is specified to sbatch
| --ntasks=1  Tells Slurm to allocate a single task
| --cpus-per-task=4  tells Slurm how many CPUs to allocate to the task
| --threads-per-core=1  specifies that one thread should be run on each physical core (i.e. ignore hyperthreading). This is usually recommended for 
cpu-intensive parallel jobs.
| --gres=lscratch:100 specifies that 100 GB of local scratch on the node should be allocated to this job. This parameter is **required** for all Psi4 job submissions,
including interactive jobs. 

 | |
 | |
 | |
 | |
 | |




The Psi4 module **must** be loaded from within a batch job. It sets the environment variable $PSI\_SCRATCH to '/lscratch/$SLURM\_JOBID'. This directory is
only created if you submit the job with --gres=lscratch:#. (see the [User Guide section on
using local disk](http://hpc.nih.gov/docs/userguide.html#local). **If you do not allocate local scratch space, your job will fail.**

Running Psi4 on the Biowulf login node or Helix is not permitted and will fail.

Documentation

[Psi4 Manual](http://psicode.org/psi4manual/master/index.html)




























