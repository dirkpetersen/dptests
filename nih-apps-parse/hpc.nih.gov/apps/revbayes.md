

document.querySelector('title').textContent = 'RevBayes on Biowulf';
RevBayes on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
[Swarm of jobs](#swarm) 
[MPI Parallel job](#MPI)
 |



 Revbayes is an application for Bayesian phylogenetic inference, and includes an MPI implementation. According to the authors:
 
> 
>  RevBayes provides an interactive environment for statistical computation in phylogenetics. It is primarily intended for modeling, simulation, and Bayesian inference in evolutionary biology, particularly phylogenetics. However, the environment is quite general and can be useful for many complex modeling tasks.
>  



> 
>  RevBayes uses its own language, Rev, which is a probabilistic programming language like JAGS, STAN, Edward, PyMC3, and related software. However, phylogenetic models require inference machinery and distributions that are unavailable in these other tools.
>  



> 
>  The Rev language is similar to the language used in R. Like the R language, Rev is designed to support interactive analysis. It supports both functional and procedural programming models, and makes a clear distinction between the two. Rev is also more strongly typed than R.
>  





### References:


* Höhna, Landis, Heath, Boussau, Lartillot, Moore, Huelsenbeck, Ronquist. 2016 
 [RevBayes: Bayesian phylogenetic inference using graphical models and an interactive model-specification language](http://sysbio.oxfordjournals.org/content/65/4/726)
*Systematic Biology*, 65:726-736.
* Höhna, Heath, Boussau, Landis, Ronquist, Huelsenbeck. 2014.
 [Probabilistic graphical model representation in phylogenetics.](http://sysbio.oxfordjournals.org/content/63/5/753)
*Systematic Biology*, 63:753–771.



Documentation
* [Revbayes main site](https://revbayes.github.io/)
* [RevBayes Tutorials](https://revbayes.github.io/tutorials/).


Important Notes
* Module Name: revbayes (see [the modules page](/apps/modules.html) for more information)
* Singlethreaded and MPI
* Environment variables set 
	+ RB\_EXAMPLE\_DATA set to /usr/local/apps/revabyes/example-data which includes tutorial data [available here](https://github.com/revbayes/revbayes_tutorial)



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

[user@cn3144 ~]$ **module load revabyes**

[user@cn3144 ~]$ **cp -r $RB\_EXAMPLE\_DATA/RB\_Partition\_Tutorial .**

[user@cn3144 ~]$ **cd RB\_Partition\_Tutorial**

[user@cn3144 ~]$ **rb scripts/mcmc\_Partition\_uniform.Rev**

RevBayes version (1.1.1)
Build from development (rapture-641-g7d4aab) on Tue Mar  9 13:42:53 EST 2021

Visit the website www.RevBayes.com for more information about RevBayes.

RevBayes is free software released under the GPL license, version 3. Type 'license()' for details.

To quit RevBayes type 'quit()' or 'q()'.


> source("scripts/mcmc_Partition_uniform.Rev")
   Processing file "scripts/mcmc_Partition_uniform.Rev"
   Successfully read one character matrix from file 'data/primates_and_galeopterus_cox2.nex'
   Successfully read one character matrix from file 'data/primates_and_galeopterus_cytb.nex'

   Running burn-in phase of Monte Carlo sampler for 10000 iterations.
   This simulation runs 2 independent replicates.
   The simulator uses 52 different moves in a random move schedule with 62 moves per iteration

Progress:
0---------------25---------------50---------------75--------------100
********************************************************************


   Running MCMC simulation
   This simulation runs 2 independent replicates.
   The simulator uses 52 different moves in a random move schedule with 62 moves per iteration

Iter        |      Posterior   |     Likelihood   | ...
----------------------------------------------------...
0           |       -21259.5   |       -21227.8   | ...
1000        |       -21273.7   |       -21241.1   | ...
[...]

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. revbayes.sh). For example:



```

#!/bin/bash
set -e
module load revbayes
cd /data/$USER/analysis_dir/
rb scripts/full_analysis.Rev

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command. Since this is the singlethreaded version, you do not need to update the number of CPUs allocated.



```
sbatch [--mem=#] revbayes.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. revbayes.swarm). For example:



```

rb script/mcmc_model1.Rev
rb script/mcmc_model2.Rev
rb script/mcmc_model3.Rev
[...]

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f revbayes.swarm [-g #] --module revbayes
```

where


|  |  |  |  |
| --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | --module revbayes Loads the revbayes module for each subjob in the swarm 
 | |
 | |



 MPI Revbayes batch job

An MPI job can run in parallel across multiple nodes. We highly recommend reading our [guidelines on multinode jobs](https://hpc.nih.gov/policies/multinode.html).
Set up a batch script, e.g. revbayes-mpi.sh, along the following lines:

```

  #!/bin/bash
  #SBATCH --partition=multinode
  #SBATCH --constraint=x2695
  #SBATCH --ntasks=8
  #SBATCH --ntasks-per-core=1
  #SBATCH --mem-per-cpu=1G
  #SBATCH --time=2:00:00
  
  module load revbayes
  
  WDIR=/data/$USER/analysis_dir
  
  cd $WDIR
  
  mpirun -np $SLURM_NTASKS rb-mpi scripts/full_analysis.Rev

```

Submit the job with:

```

sbatch revbayes-mpi.sh

```













