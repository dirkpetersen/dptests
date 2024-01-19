

document.querySelector('title').textContent = 'APBS (Adaptive Poisson-Boltzmann Solver) on Biowulf';
APBS (Adaptive Poisson-Boltzmann Solver) on Biowulf
APBS is a software package for the numerical solution of the
Poisson-Boltzmann equation (PBE), one of the most popular continuum models for
describing electrostatic interactions between molecular solutes in salty,
aqueous media. Continuum electrostatics plays an important role in several
areas of biomolecular simulation, including:


* simulation of diffusional processes to determine ligand-protein and
protein-protein binding kinetics,
* implicit solvent molecular dynamics of biomolecules ,
* solvation and binding energy calculations to determine ligand-protein and
protein-protein equilibrium binding constants and aid in rational drug
design,
* and biomolecular titration studies.


APBS was designed to efficiently evaluate electrostatic properties for such
simulations for a wide range of length scales to enable the investigation of
molecules with tens to millions of atoms.


There are multiple versions of APBS available. An easy way of selecting the version is to use [modules](/apps/modules.html). To see the modules available, type



```
module avail apbs
```

To select a module, type



```
module load apbs/[ver]
```

where [ver] is the version of choice. This will set your $PATH variable to allow the apbs executables.


Please note that the APBS 3.0.0 module contains pdb2pqr version 3.4.0.


Interactive use
---------------



```
module load apbs
apbs < apbs.in > apbs.out
```

sbatch
------


Create a batch input file (e.g. apbs\_run.sh), which uses the input file
'apbs.in'. For example:



```
#!/bin/bash
module load apbs
apbs < apbs.in > apbs.out
```

Submit this job using the Slurm 'sbatch' command.



```
sbatch --cpus-per-task=1 apbs-run.sh
```

Examples
--------


Input files are available in $APBSEXAMPLES:



```
ls $APBSEXAMPLES
actin-dimer  bem   FKBP     hca-bind  ion-pmf      membrane  opal     point-pmf    README.html  solv
alkanes      born  geoflow  ionize    ion-protein  misc      pka-lig  protein-rna  smpbe

```

Documentation
-------------


* Main web site: <http://www.poissonboltzmann.org/>
* APBS Users Mailing List: <http://lists.sourceforge.net/mailman/listinfo/apbs-users>
- an excellent resource for detailed discussions of the software




