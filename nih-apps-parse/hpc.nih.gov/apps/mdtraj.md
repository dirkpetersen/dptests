

document.querySelector('title').textContent = 'mdtraj on Biowulf';
mdtraj on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
[Swarm of jobs](#swarm) 
 |



MDTraj is a python library that allows users to manipulate molecular dynamics (MD) trajectories and perform a variety of analyses, including fast RMSD, solvent accessible surface area, hydrogen bonding, etc. A highlight of MDTraj is the wide variety of molecular dynamics trajectory file formats which are supported, including RCSB pdb, GROMACS xtc and trr, CHARMM / NAMD dcd, AMBER binpos, AMBER NetCDF, AMBER mdcrd, TINKER arc and MDTraj HDF5. 




The installation of MDTraj on Helix & Biowulf also contains MDAnalysis, another python package for analyzing MD trajectories. 



### References:


* McGibbon RT, Beauchamp KA, Harrigan MP, Klein C, Swails JM, Hernández CX, Schwantes CR, Wang LP, Lane TJ, Pande VS.
 [**MDTraj: A Modern Open Library for the Analysis of Molecular Dynamics Trajectories.**](https://www.ncbi.nlm.nih.gov/pubmed/26488642)
*Biophys J. 2015 Oct 20;109(8):1528-32.*


Documentation
* [MDTraj Documentation](http://mdtraj.org)
* [MDAnalysis Documentation](http://www.mdanalysis.org//)


Important Notes
* Module Name: mdtraj (see [the modules page](/apps/modules.html) for more information)
* Singlethreaded
* Environment variables set 
	+ MDTRAJ\_EXAMPLES* Example files in $MDTRAJ\_EXAMPLES



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

[user@cn3144 ~]$ module load mdtraj
[user@cn3144 ~]$ python
>>> import mdtraj as md
>>> t = md.load('trajectory.xtc', top='trajectory.pdb')
>>> print t
<mdtraj.Trajectory with 100 frames, 22 atoms at 0x109f0a3d0>
>>> quit()
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a python input script **mdtraj.py**. This is an example lifted from the
[MDTraj Introduction](http://mdtraj.org/latest/examples/introduction.html) page.



```

#!/usr/bin/env python
import mdtraj as md
t = md.load('trajectory.xtc', top='trajectory.pdb')
print t
print t[1:10]
print t[-1]
print t.xyz.shape
print np.mean(t.xyz)
print t.time[0:10]
t.unitcell_lengths[-1]
t[::2].save('halftraj.h5')
t[0:10].save_dcd('first-ten-frames.dcd')
atom_to_keep = [a.index for a in t.topology.atoms if a.name == 'CA']
t.restrict_atoms(atoms_to_keep)  # this acts inplace on the trajectory
t.save('CA-only.h5')

```

Create a batch script **mdtraj.sh**:



```
#!/bin/bash
module load mdtraj
python mdtraj.py
```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] mdtraj.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. mdtraj.swarm). For example:



```

python mdtraj_1.py
python mdtraj_2.py
python mdtraj_3.py
python mdtraj_4.py

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f mdtraj.swarm [-g #] [-t #] --module mdtraj
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module mdtraj Loads the mdtraj module for each subjob in the swarm 
 | |
 | |
 | |








