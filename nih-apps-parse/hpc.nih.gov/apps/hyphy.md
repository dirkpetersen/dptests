

document.querySelector('title').textContent = 'HYPHY on Biowulf';
HYPHY on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job (threaded version)](#int-threaded)
[Interactive job (MPI version)](#int-mpi)
[Batch job -- threaded](#sbatch-threaded)
[Batch job -- MPI](#sbatch-MPI)
 |



HyPhy (Hypothesis Testing using Phylogenies) is an open-source software package for the analysis of genetic sequences (in particular the inference of natural selection) using techniques in phylogenetics, molecular evolution, and machine learning. It features a rich scripting language for limitless customization of analyses. 



Documentation
* [HYPHY Main Site](http://hyphy.org/)
* [Command-line tool documentation](http://hyphy.org/tutorials/CLI-tutorial/)


Important Notes
* Module Name: hyphy (see [the modules page](/apps/modules.html) for more information)
* Multithreaded/MPI
* Environment variables set
+ HYPHY\_TEMPLATES -- the location for the hyphy template .bf files
+ HYPHY\_TESTDATA -- the location for the test data used in the examples below


HYPHYMPI v2.5.29 does not work on e7543 nodes. See the [MPI Interactive session](#int-mpi) below for an example of how to allocate a different node type



Interactive job 
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive --cpus-per-task=8 --mem=20g --gres=lscratch:5**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load hyphy**

[user@cn3144 ~]$ **cd /lscratch/$SLURM\_JOB\_ID**

[user@cn3144 ~]$ **cp $HYPHY\_TESTDATA/CD2.\* .**

[user@cn3144 ~]$ **hyphy CPU=$SLURM\_CPUS\_PER\_TASK GARD --alignment CD2.nex --tree CD2.newick**

Analysis Description
--------------------
GARD : Genetic Algorithms for Recombination Detection. Implements a
heuristic approach to screening alignments of sequences for
recombination, by using the CHC genetic algorithm to search for
phylogenetic incongruence among different partitions of the data. The
number of partitions is determined using a step-up procedure, while the
placement of breakpoints is searched for with the GA. The best fitting
model (based on c-AIC) is returned; and additional post-hoc tests run to
distinguish topological incongruence from rate-variation. v0.2 adds and
spooling results to JSON after each breakpoint search conclusion

- __Requirements__: A sequence alignment.

- __Citation__: **Automated Phylogenetic Detection of Recombination Using a Genetic
Algorithm**, _Mol Biol Evol 23(10), 1891–1901

- __Written by__: Sergei L Kosakovsky Pond

- __Contact Information__: spond@temple.edu

- __Analysis Version__: 0.2

type: nucleotide
rv: None
>Loaded a nucleotide multiple sequence alignment with **10** sequences, **561** sites (390 of which are variable) from `/lscratch/46116226/CD2.nex`
>Minimum size of a partition is set to be 17 sites


### Fitting the baseline (single-partition; no breakpoints) model
* Log(L) = -3529.89, AIC-c =  7112.21 (25 estimated parameters)

### Performing an exhaustive single breakpoint analysis
Done with single breakpoint analysis.
   Best sinlge break point location: 25
   c-AIC  = 7106.656403051507

### Performing multi breakpoint analysis using a genetic algorithm
Done with 2 breakpoint analysis.
    Best break point locations: 25, 65
    c-AIC = 7088.977422652543
Done with 3 breakpoint analysis.
    Best break point locations: 25, 65, 173
    c-AIC = 7101.784222144825

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```

Interactive job - MPI version
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) with ntasks and run the program with mpirun.


HYPHYMPI does not run on e7543 nodes. Use the [--constraint](https://slurm.schedmd.com/sbatch.html#OPT_constraint) option to choose another node type like x6140 or x2695. See the feature table in the output of freen command to find different node types. 


Sample session (user input in **bold**):




```

[user@biowulf]$ **sinteractive --ntasks=8 --ntasks-per-core=1 --mem=20g --constraint=x6140**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load hyphy**

[user@cn3144 ~]$ **mkdir /data/$USER/hyphy**

[user@cn3144 ~]$ **cp $HYPHY\_TESTDATA/CD2.\* /data/$USER/hyphy/** 

[user@cn3144 ~]$ **mpirun -np $SLURM\_NTASKS HYPHYMPI gard** 


Analysis Description
--------------------
GARD : Genetic Algorithms for Recombination Detection. Implements a
heuristic approach to screening alignments of sequences for
recombination, by using the CHC genetic algorithm to search for
phylogenetic incongruence among different partitions of the data. The
number of partitions is determined using a step-up procedure, while the
placement of breakpoints is searched for with the GA. The best fitting
model (based on c-AIC) is returned; and additional post-hoc tests run to
distinguish topological incongruence from rate-variation.

- __Requirements__: A sequence alignment.

- __Citation__: **Automated Phylogenetic Detection of Recombination Using a Genetic
Algorithm**, _Mol Biol Evol 23(10), 1891-1901

- __Written by__: Sergei L Kosakovsky Pond

- __Contact Information__: spond@temple.edu

- __Analysis Version__: 0.1

type: Nucleotide

Select a sequence alignment file **/data/$USER/hyphy/CD2.nex**
rv: None
>Loaded a Nucleotide multiple sequence alignment with **10** sequences, **561** sites (390 of which are variable) from `/data/$USER/hyphy/CD2.nex`
>Minimum size of a partition is set to be 17 sites


### Fitting the baseline (single-partition; no breakpoints) model
* Log(L) = -3529.89, AIC-c =  7112.21 (25 estimated parameters)

### Performing an exhaustive single breakpoint analysis
Done with single breakpoint analysis.
   Best sinlge break point location: 25
   c-AIC  = 7107.094096738777

### Performing multi breakpoint analysis using a genetic algorithm
[...]

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$


```


Batch job - threaded version
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. hyphy.sh). For example, to run the threaded version of hyphy:



```

#!/bin/bash
set -e
module load hyphy
cp $HYPHY_TESTDATA/CD2* .

hyphy CPU=$SLURM_CPUS_PER_TASK slac --alignment CD2.nex --tree CD2.newick

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=# [--mem=#] hyphy.sh
```


Batch job - MPI version

To run the MPI version of hyphy, here is a sample batch script. 

```

#!/bin/bash

module load hyphy
# copy the test data
mkdir /data/$USER/hyphy
cp $HYPHY_TESTDATA/CD2.*   /data/$USER/hyphy/ 

mpirun -np $SLURM_NTASKS HYPHYMPI gard --alignment CD2.nex

```

Submit this job with:

```

sbatch --ntasks=8 --ntasks-per-core=1 --mem=20g  --constraint=x6140 hyphympi.bat

```













