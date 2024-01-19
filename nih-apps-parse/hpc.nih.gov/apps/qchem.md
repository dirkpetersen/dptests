

document.querySelector('title').textContent = 'Q-Chem on Biowulf';
Q-Chem on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[OpenMP (threaded) Batch job](#openmp) 
[MPI batch job](#mpi)
 |



Q-Chem is a comprehensive ab initio quantum chemistry package for accurate predictions of molecular structures, reactivities, and vibrational, electronic and NMR spectra. The new release of Q-Chem 4.0 represents the state-of-the-art of methodology from the highest performance DFT/HF calculations to high level post-HF correlation methods:

* Dispersion-corrected and double hybrid DFT functionals
* Faster algorithms for DFT, HF and coupled-cluster calculations;
* Structures and vibrations of excited states with TD-DFT;
* Methods for mapping complicated potential energy surfaces;
* Efficient valence space models for strong correlation;
* More choices for excited states, solvation and charge-transfer;
* Effective Fragment Potential and QM/MM for large systems;
* Shared-memory for multicores and implementations for GPU's.





Q-Chem is a licensed product developed by [Q-Chem](http://www.q-chem.com)
Documentation
* [Q-Chem 5.0 manual](http://www.q-chem.com/qchem-website/manual/qchem50_manual/index.html)


Important Notes
* Module Name: qchem (see [the modules page](/apps/modules.html) for more information)
* Q-Chem can use threads or MPI processes to parallelize. Different methods within Q-Chem can utilize one or both kinds of parallelization. See the
section on [Running Q-Chem in parallel](http://www.q-chem.com/qchem-website/manual/qchem50_manual/sect-running.html) for
more information. 
**OpenMP (threads)** can utilize multiple cores on a node, but cannot run on multiple nodes. 
* After you run 'module load qchem', a set of example input files is available in $QC/samples.
* Reference data in /fdb/qchem/



OpenMP (theaded) Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. qchem.sh). For example:



```

#!/bin/bash
#SBATCH --job-name="QC"
#SBATCH --mail-type=BEGIN,END

module load qchem/4.3

cd /data/$USER/qchem

# copy the sample data to this directory
cp ${QC}/samples/bsse/frgm_mp2_h2o_h2o_h2o.in .

# run qchem
qchem -nt $SLURM_CPUS_PER_TASK  frgm_mp2_h2o_h2o_h2o.in   frgm_mp2_h2o_h2o_h2o.in.out

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command. e.g.



```
sbatch --cpus-per-task=16  --threads-per-core=1 --mem=20g qchem.sh
```



|  |  |  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- | --- | --- |
| -nt $SLURM\_CPUS\_PER\_TASK  tells Q-Chem to run $SLURM\_CPUS\_PER\_TASK threads
| --cpus-per-task=16  tells Slurm how many CPUs to allocate
| --threads-per-core=1  specifies that one thread should be run on each physical core (i.e. ignore hyperthreading). This is usually recommended for 
cpu-intensive parallel jobs.
| --gres=lscratch:100 specifies that 100 GB of local scratch on the node should be allocated to this job. This parameter is required for all Qchem job submissions. 
 | |
 | |
 | |
 | |



When the qchem module is loaded as part of a batch job, it sets the environment variable $QCLOCALSCR to '/lscratch/$SLURM\_JOBID' . This directory is
only created if you submit the job with --gres=lscratch:#. (see the [User Guide section on
using local disk](http://hpc.nih.gov/docs/userguide.html#local). 

MPI batch job on Biowulf 

Q-Chem can also run **MPI** . Sample MPI batch script

```


#!/bin/bash
#SBATCH --job-name="QC"

module load qchem/4.3

cd /data/$USER/qchem

# copy the sample input file to this directory
cp  ${QC}/samples/sp/dft_b5050lyp_c4h6.in .

# run qchem
make-qchem-nodefile
export QCMACHINEFILE=`pwd`/qchem.$SLURM_JOBID.nodes
qchem -np $SLURM_NTASKS  dft_b5050lyp_c4h6.in   dft_b5050lyp_c4h6.in.out
rm $QCMACHINEFILE


```


Submit this job with:

```

biowulf% sbatch --ntasks=4 --ntasks-per-core=1  --gres=lscratch:100 qchem.bat

```


where:


|  |  |  |  |  |  |  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| make-qchem-nodefile  a script that writes a file containing a list of nodes to be used
| $QCMACHINEFILE  environment variable used by Q-Chem to determine the nodes on which to run
| -np $SLURM\_NTASKS  tells Q-Chem to run $SLURM\_NTASKS MPI processes
| --ntasks=4  tells Slurm how many tasks (MPI processes) to run
| --ntasks-per-core=1  tells Slurm to run only one task on each physical core (i.e. ignore hyperthreading). This is usually recommended for 
cpu-intensive parallel jobs.
| --gres=lscratch:100 specifies that 100 GB of local scratch on the node should be allocated to this job. This parameter is required for all Qchem job submissions. 
 | |
 | |
 | |
 | |
 | |
 | |












