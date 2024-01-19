

document.querySelector('title').textContent = 'Rosetta on Biowulf';
Rosetta on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
[Multithreaded job](#multithreaded) 
[MPI job](#mpi) 
[Swarm of jobs](#swarm) 
 |



The Rosetta++ software suite focuses on the prediction and design of protein structures, protein folding mechanisms, and protein-protein interactions. The Rosetta codes have been repeatedly successful in the Critical Assessment of Techniques for Protein Structure Prediction (CASP) competition as well as the CAPRI competition and have been modified to address additional aspects of protein design, docking and structure.



Documentation
* Rosetta Main Site: <http://www.rosettacommons.org/>


Important Notes
* Module Name: rosetta (see [the modules page](/apps/modules.html) for more information)
* Multithreaded/singlethreaded/MPI
* Unusual environment variables set:
	+ **ROSETTA3\_HOME** -- path to the Rosetta installation directory
	+ **ROSETTA3\_DB** -- path to the Rosetta database directory
	+ **ROSETTA\_DATABASE** -- ditto
	+ **ROSETTA\_TOOLS** -- path to miscellaneous scripts and wrappers* Example script in $ROSETTA3\_HOME/../rosetta3\_demos.tgz


**NOTE FOR LARGE ROSETTA JOBS:** When running jobs > 50 cpus, it is best practice to copy the Rosetta database to local scratch, using the following command:



```

sbcast ${ROSETTA3_DB}.tgz /lscratch/$SLURM_JOB_ID/database.tgz && srun tar -C /lscratch/$SLURM_JOB_ID/ -xzf /lscratch/$SLURM_JOB_ID/database.tgz
export ROSETTA3_DB=/lscratch/$SLURM_JOB_ID/database
export ROSETTA_DATABASE=/lscratch/$SLURM_JOB_ID/database

```

If the -database is used as an option, you must specify the location as -database /lscratch/$SLURM\_JOB\_ID/database. In addition, you must allocate at least 10 GB of local scratch space with --gres=lscratch:10.



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

[user@cn3144 ~]$ module load rosetta
[user@cn3144 ~]$ tar xzvf $ROSETTA3_HOME/../rosetta3_demos.tgz
[user@cn3144 ~]$ ./run_demo.sh

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```

Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. rosetta.sh). For example:



```

#!/bin/bash
module load rosetta
relax @flags > relax.log

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--mem=#] rosetta.sh
```


Multithreaded job
A small fraction of Rosetta commands can be accelerated using multithreading. This requires an extra option:



```
-multithreading:total_threads ***#***
```

where ***#*** is the number of threads.


The easiest way to handle this is to supply a slurm-created environment variable **$SLURM\_CPUS\_ON\_NODE** in a job with multiple cpus allocated. Here is an example sbatch submission script:



```
#!/bin/bash
module load rosetta
relax @flags -multithreading:total_threads $SLURM_CPUS_ON_NODE > relax.log

```

The multiple cpus are allocated like so:



```
[user@biowulf ~]$ sbatch **--cpus-per-task=#** rosetta.sh
```

where # is a number between 2 and 8.


**NOTE:** There are very few rosetta commands that benefit from multithreading. Make sure your jobs actually utilize the excess cpus using [**jobload**](https://hpc.nih.gov/docs/biowulf_tools.html#jobload) or the [**HPC dashboard**](https://hpc.nih.gov/dashboard).



MPI job
Certain Rosetta commands can be distributed using MPI. This requires the MPI module for Rosetta be loaded, e.g.



```
module load rosetta/2018.21.mpi
```

The Rosetta commands are launched using **mpirun**, passing the number of tasks with **-np**:



```
mpirun -np ${SLURM_NTASKS} AbinitioRelax @flags
```

Then, tasks must be allocated instead of CPUs:



```
sbatch [--ntasks=#] [--ntasks-per-core=1] rosetta.sh
```

Not all Rosetta commands are MPI-enabled. Check the [Rosetta Commons documentation](https://www.rosettacommons.org/docs/latest/rosetta_basics/MPI) to learn more about running Rosetta with MPI.


Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. rosetta.swarm). For example:



```
AbinitioRelax @flags -out:file:silent abinito1.out > abinitio1.log
AbinitioRelax @flags -out:file:silent abinito2.out > abinitio2.log
AbinitioRelax @flags -out:file:silent abinito3.out > abinitio3.log
AbinitioRelax @flags -out:file:silent abinito4.out > abinitio4.log
AbinitioRelax @flags -out:file:silent abinito5.out > abinitio5.log
AbinitioRelax @flags -out:file:silent abinito6.out > abinitio6.log
AbinitioRelax @flags -out:file:silent abinito7.out > abinitio7.log
AbinitioRelax @flags -out:file:silent abinito8.out > abinitio8.log

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f rosetta.swarm [-g #] [-t #] --module rosetta
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g #  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t #  Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module rosetta  Loads the rosetta module for each subjob in the swarm 
 | |
 | |
 | |






