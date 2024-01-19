

document.querySelector('title').textContent = 'Phenix on Biowulf';
Phenix on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
[Swarm of jobs](#swarm) 
[Running Rosetta through Phenix](#rosetta) 
[CryoFit](#cryofit)
 |



![Phenix logo](/images/phenix.jpg)
The [PHENIX software suite](http://www.phenix-online.org) is a highly automated system for macromolecular structure determination that can rapidly arrive at an initial partial model of a structure without significant human intervention, given moderate resolution and good quality data.





### References:


* P.D. Adams, P.V. Afonine, G. Bunkoczi, V.B. Chen, I.W. Davis, N. Echols, J.J. Headd, L.W. Hung, G.J. Kapral, R.W. Grosse-Kunstleve, A.J. McCoy, N.W. Moriarty, R. Oeffner, R.J. Read, D.C. Richardson, J.S. Richardson, T.C. Terwilliger, and P.H. Zwart.
 [**PHENIX: a comprehensive Python-based system for macromolecular structure solution.**](http://www.ncbi.nlm.nih.gov/pubmed/21821126)
*Acta Cryst. D66, 213-221 (2010)*


Documentation
* [Phenix Documentation](http://www.phenix-online.org/documentation/)
* [Crystallographic Refinement with Rosetta](https://www.phenix-online.org/documentation/reference/rosetta_refine.html)
* [Fit Biomolecules into Cryo-EM Maps using MD Simulation](https://phenix-online.org/documentation/tutorials/cryo_fit_gui.html)


Important Notes
* Module Name: Phenix (see [the modules page](/apps/modules.html) for more information)
* Singlethreaded
* Environment variables set 
	+ **$PATH**
	+ **$PHENIX** -- main installation directory
	+ **$PHENIX\_VERSION** -- version
	+ **$PHENIX\_MODULES\_DIR** -- path to Phenix modules
	+ **$PHENIX\_ROSETTA\_PATH** -- Rosetta installation directory
	+ **$ROSETTA3\_DB** -- Rosetta database directory* Extra modules:
	+ **[cryo\_fit](https://phenix-online.org/documentation/reference/cryo_fit.html)** -- fit biomolecule to a cryoem map


This application requires an [X-Windows connection](/docs/connect.html). Further, there are combinations of X11 servers and drivers that cause Phenix to crash. It is known that XQuartz (v2.7.x) is incompatible with Phenix. Users are encouraged to use [NX or FastX](https://hpc.nih.gov/docs/nx.html) as their X11 servers.


Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load Phenix**
[user@cn3144 ~]$ **phenix**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```

### Pymol


Phenix can launch graphical applications Coot and Pymol. Navigate to Preferences -> Graphics and insert this string into the text box for PyMOL path:


* **Coot:** /usr/local/apps/Coot/0.9.8.92/bin/coot
* **Pymol:** /usr/local/apps/pymol/2.1.0/bin/pymol


![PyMOL path](Phenix.jpg)
Now when the PyMOL button is clicked, you should see this:


![PyMOL window](PyMOL.png)

Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. Phenix.sh). For example:



```

#!/bin/bash
module load Phenix
phenix.autobuild seq_file=p9.seq data=p9-solve.mtz \
  input_map_file=p9-resolve.mtz resolution=2.4  \
  ncs_copies=1 nproc=$SLURM_CPUS_ON_NODE \
  temp_dir=/lscratch/$SLURM_JOB_ID
```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] --gres=lscratch:50 Phenix.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. Phenix.swarm). For example:



```

phenix.elbow input_1.pdb
phenix.elbow input_2.pdb
phenix.elbow input_3.pdb
phenix.elbow input_4.pdb
```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f Phenix.swarm [-g #] [-t #] --module Phenix
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module Phenix Loads the Phenix module for each subjob in the swarm 
 | |
 | |
 | |


Using Rosetta within Phenix
[Rosetta](rosetta.html) has been compiled specifically for Phenix, and is available via
two environment variables **PHENIX\_ROSETTA\_PATH** and **ROSETTA3\_DB**. Do NOT load the Rosetta modules built for
Biowulf, as these are NOT compiled specifically for Phenix. Also note that structural models are limited to standard amino acids and other limitations.


Here is an example of structural refinement using Rosetta to generate structural models (this file is named Rosetta\_refine.sh):



```
#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --time=72:00:00
#SBATCH --mem=50g
#SBATCH --gres=lscratch:50
module load Phenix
phenix.rosetta_refine input.pdb input.mtz \
  nproc=$SLURM_CPUS_ON_NODE \
  temp_dir=/lscratch/$SLURM_JOB_ID \
  post_refine=True
```

Then submit to the cluster:



```
sbatch Rosetta_refine.sh
```

Using CryoFit within Phenix
The [CryoFit](https://phenix-online.org/documentation/tutorials/cryo_fit_gui.html) protocol can be used to fit a fasta or pdb into a cryoEM map.


* **gromacs\_cryo\_fit:** /usr/local/apps/Phenix/1.20.1-4487/cryo\_fit-2Jan2020/cryo\_fit/bin


Both cryo\_fit and cryo\_fit2 can be run from the commandline. Here is how to run the tutorials:



```

cp $PHENIX_MODULES_DIR/cryo_fit/tutorial_input_files/*.{mrc,pdb} .
phenix.cryo_fit GTPase_activation_center_tutorial.pdb GTPase_activation_center_tutorial_gaussian_1p5.mrc

```


```

cp -R $PHENIX_MODULES_DIR/cryo_fit2/tutorial/input .
phenix.cryo_fit2 input/tutorial_cryo_fit2_model.pdb input/tutorial_cryo_fit2_map.ccp4 resolution=4 

```





