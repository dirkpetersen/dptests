

document.querySelector('title').textContent = 'Open Babel on Biowulf';
Open Babel on Biowulf
![Open Babel Logo](OpenBabel.png)


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


Open Babel is a chemical toolbox designed to speak the many languages of chemical data. It's an open, collaborative project allowing anyone to search, convert, analyze, or store data from molecular modeling, chemistry, solid-state materials, biochemistry, or related areas. 


### References:


* Noel M O'Boyle, Michael Banck, Craig A James, Chris Morley, Tim Vandermeersch and Geoffrey R Hutchison
 [**Open Babel: An open chemical toolbox**](https://jcheminf.springeropen.com/articles/10.1186/1758-2946-3-33)
*Journal of Cheminformatics (2011) 3:33*


Documentation
* [Open Babel Wiki Site](http://openbabel.org)


Important Notes
* Module Name: OpenBabel (see [the modules page](/apps/modules.html) for more information)
 * Singlethreaded
 * Environment variables set 
	+ OPENBABEL\_HOME
	+ OPENBABEL\_EXAMPLES* Example files in $OPENBABEL\_EXAMPLES


The executables are version dependent:



* 2.4.1
* 3.0.0



### Executables


* **[babel](http://openbabel.org/wiki/Babel)** - a converter for chemistry and molecular modeling data files
* **[obabel](https://openbabel.org/docs/dev/Command-line_tools/babel.html)** - a converter for chemistry and molecular modeling data files
* **[obchiral](http://openbabel.org/wiki/Obchiral)** - print molecular chirality information
* **[obconformer](http://openbabel.org/wiki/Obconformer)** - generate low-energy conformers
* **[obenergy](http://openbabel.org/wiki/Obenergy)** - calculate the energy for a molecule
* **[obfit](http://openbabel.org/wiki/Obfit)** - superimpose two molecules based on a SMARTS pattern
* **[obgen](http://openbabel.org/wiki/Obgen)** - generate 3D coordinates for a molecule
* **[obgrep](http://openbabel.org/wiki/Obgrep)** - an advanced molecular grep program using SMARTS
* **[obminimize](http://openbabel.org/wiki/Obminimize)** - optimize the geometry, minimize the energy for a molecule
* **[obprobe](http://openbabel.org/wiki/Obprobe)** - create electrostatic probe grid
* **[obprop](http://openbabel.org/wiki/Obprop)** - print standard molecular properties
* **[obrotamer](http://openbabel.org/wiki/Obrotamer)** - generate conformer/rotamer coordinates
* **[obrotate](http://openbabel.org/wiki/Obrotate)** - batch-rotate dihedral angles matching SMARTS patterns
* **[obspectrophore](http://openbabel.org/docs/current/Fingerprints/spectrophore.html)** - calculate a Spectrophore, a one-dimensional descriptor generated from the property field surrounding the molecule


 

### Executables


* **[obabel](https://openbabel.org/docs/dev/Command-line_tools/babel.html)** - a converter for chemistry and molecular modeling data files
* **[obconformer](http://openbabel.org/wiki/Obconformer)** - generate low-energy conformers
* **[obenergy](http://openbabel.org/wiki/Obenergy)** - calculate the energy for a molecule
* **[obfit](http://openbabel.org/wiki/Obfit)** - superimpose two molecules based on a SMARTS pattern
* **obfitall** - ???
* **[obgen](http://openbabel.org/wiki/Obgen)** - generate 3D coordinates for a molecule
* **[obgrep](http://openbabel.org/wiki/Obgrep)** - an advanced molecular grep program using SMARTS
* **[obminimize](http://openbabel.org/wiki/Obminimize)** - optimize the geometry, minimize the energy for a molecule
* **obmm** - OpenBabel molecular mechanics program
* **[obprobe](http://openbabel.org/wiki/Obprobe)** - create electrostatic probe grid
* **[obprop](http://openbabel.org/wiki/Obprop)** - print standard molecular properties
* **[obrotamer](http://openbabel.org/wiki/Obrotamer)** - generate conformer/rotamer coordinates
* **[obrotate](http://openbabel.org/wiki/Obrotate)** - batch-rotate dihedral angles matching SMARTS patterns
* **[obspectrophore](http://openbabel.org/docs/current/Fingerprints/spectrophore.html)** - calculate a Spectrophore, a one-dimensional descriptor generated from the property field surrounding the molecule


 

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

[user@cn3144 ~]$ **module load OpenBabel**
[user@cn3144 ~]$ **obabel -ixyz benzene.xyz -O benzene.pdb**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. OpenBabel.sh). For example:



```

#!/bin/bash
module load OpenBabel
obabel -ipdb benzene.pdb -O benzene.xyz

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] OpenBabel.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. OpenBabel.swarm). For example:



```

obabel mymols.sdf -f 1 -l 10 -Ooutputfile1.sdf
obabel mymols.sdf -f 11 -l 20 -O outputfile2.sdf
obabel mymols.sdf -f 21 -l 30 -O outputfile3.sdf
obabel mymols.sdf -f 31 -l 40 -O outputfile4.sdf

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f OpenBabel.swarm --module OpenBabel
```

 where
 

|  |  |
| --- | --- |
| --module OpenBabel Loads the OpenBabel module for each subjob in the swarm 
  | |






