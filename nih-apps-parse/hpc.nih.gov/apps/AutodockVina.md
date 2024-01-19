

document.querySelector('title').textContent = 'Autodock & AutodockVina on Biowulf';
Autodock & AutodockVina on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
[Swarm of jobs](#swarm) 
[Autodock-GPU](#gpu)
 |



AutoDock is a suite of automated docking tools. It is designed to predict how small molecules, such as substrates or drug candidates, bind to a receptor of known 3D structure.

AutoDock Vina does not require choosing atom types and pre-calculating grid maps for them. Instead, it calculates the grids internally, for the atom types that are needed, and it does this virtually instantly.

Documentation
* [Autodock website at Scripps](http://autodock.scripps.edu/)


Important Notes
* Module Name: Autodock or AutodockVina (see [the modules page](/apps/modules.html) for more information)
* AutodockVina is multithreaded. Autodock is single-threaded.
* The mgltools package that is used to create input files for Vina is also installed. When you load the Autodock or AutodockVina modules, you will also get access to the
 mgltools executables.
* If you wish to run the mgltools Python scripts you may need to use the `pythonsh` command available when you load the module. This will ensure your environment is properly set up.
* mgltools is a graphics-intensive package, and can easily be installed/run on your desktop system for better performance. It can be [downloaded here](http://mgltools.scripps.edu/downloads).
* The entire PDB (updated once a week) database is available in /pdb on all nodes, the Biowulf login node, and Helix.



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --cpus-per-task=8**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load AutodockVina**

[user@cn3144 ~]$ **module list**

Currently Loaded Modules:
  1) AutodockVina/1.1.2

[user@cn3144 ~]$ **vina --config conf.txt --cpus $SLURM\_CPUS\_PER\_TASK**

#################################################################
# If you used AutoDock Vina in your work, please cite:          #
#                                                               #
# O. Trott, A. J. Olson,                                        #
# AutoDock Vina: improving the speed and accuracy of docking    #
# with a new scoring function, efficient optimization and       #
# multithreading, Journal of Computational Chemistry 31 (2010)  #
# 455-461                                                       #
#                                                               #
# DOI 10.1002/jcc.21334                                         #
#                                                               #
# Please see http://vina.scripps.edu for more information.      #
#################################################################

Output will be ligand_out.pdbqt
Reading input ... done.
Setting up the scoring function ... done.
Analyzing the binding site ... done.
Using random seed: 653328095
Performing search ...
0%   10   20   30   40   50   60   70   80   90   100%
|----|----|----|----|----|----|----|----|----|----|
***************************************************
done.
Refining results ... done.

mode |   affinity | dist from best mode
     | (kcal/mol) | rmsd l.b.| rmsd u.b.
-----+------------+----------+----------
   1        -12.3      0.000      0.000
   2        -10.1      5.153      9.937
   3         -9.6      5.800      9.437
   4         -9.5      5.027     10.364
Writing output ... done.

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. autodock.sh). For example:



```

#!/bin/bash
module load Autodock
autodock4 -p myfile.dpf -l myfile.dlg 

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--mem=#] autodock.sh
```

The --mem=#g flag should be used if the autodock run requires more than the default 4 GB of memory.

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. myfile.swarm). For example:



```

cd /data/$USER/mydir;  autodock4 -p lig1.macro.def -l lig1.log
cd /data/$USER/mydir;  autodock4 -p lig2.macro.def -l lig2.log
cd /data/$USER/mydir;  autodock4 -p lig3.macro.def -l lig3.log
cd /data/$USER/mydir;  autodock4 -p lig4.macro.def -l lig4.log
[...]

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f myfile.swarm [-g #] --module Autodock
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g #  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t #  Number of threads/CPUs required for each process (1 line in the swarm command file). (not useful for Autodock which is single-threaded, 
 but can be used for a swarm of Vina jobs)
 | --module Autodock Loads the Autodock module for each subjob in the swarm 
 | |
 | |
 | |



Autodock - GPU

Autodock-GPU processes ligand-receptor poses in parallel over multiple compute units on GPUs. [[Autodock-GPU website](https://github.com/ccsb-scripps/AutoDock-GPU)]

Test script: 

```

#!/bin/bash

module load Autodock-GPU
cd /data/$USER
mkdir autodock-gpu
cd autodock-gpu
git clone https://github.com/L30nardoSV/reproduce-parcosi-moleculardocking.git
cd reproduce-parcosi-moleculardocking
git clone https://gitlab.com/L30nardoSV/ad-gpu_miniset_20.git

# note: this test job requires the executables to be available in the same dir
ln -s ${AUTODOCK_BIN}/* .

./prepare_inputs.sh 
./evaluate_numwi.sh <<EOF
Y
Y
1
k80
EOF

```


Select a GPU type (e.g. k80, p100, v100) and submit with:

```

sbatch --partition=gpu --gres=gpu:k80:1  test.sh

```

















