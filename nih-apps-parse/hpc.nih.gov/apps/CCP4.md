

document.querySelector('title').textContent = 'CCP4 on Biowulf';
CCP4 on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
 |


CCP4 (The Collaborative Computational Project Number 4 in Protein
 Crystallography) is a suite of separate programs which communicate via
 standard data files, rather than all operations being integrated into one
 huge program. The programs range from coordinate file manipulations and
 dataset management to sequence alignment and structure validation.


To use the programs the user must assign input and output files, and
 run the programs. Often an output file becomes the input to the next step,
 and system parameter substitution may be used to create filenames in a
 systematic way. Most crystallographic calculations involve a series of
 steps in which no decisions need be made until the end, and a command
 file provides an easy way of chaining calculations.


Standard file formats are defined for the principal data types used in
 crystallography: reflection data; density maps; and atom coordinates.
 In defining these formats, a number of trade-offs are made between
 efficiency (in space and access time), flexibility, portability, and
 simplicity of use.


There is a policy of continual technical and scientific updates to the
 suite. Where existing programs have been incorporated into the suite they
 have often subsequently undergone considerable modification above that
 needed to use the CCP4 file formats.


### References:


* Please see <http://www.ccp4.ac.uk/dist/html/REFERENCES.html> for details on citing CCP4 programs.


Documentation
* [CCP4 Documentation](http://www.ccp4.ac.uk/docs.php)


Important Notes
* Module Name: CCP4 (see [the modules page](/apps/modules.html) for more information)
* Multithreaded/Singlethreaded
* environment variables set 
	+ CCP4\_SCR -- scratch directory for CCP4, set to /lscratch/$SLURM\_JOB\_ID for interactive and batch jobs.* Reference data in /pdb/


CCP4 requires the allocation of local scratch space. Please see
<https://hpc.nih.gov/docs/userguide.html#local>
for more information.



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
This example aligns 3ESZ onto 3F6R, based on the superposition of the
 A chains from both PDB entries.


Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --gres=lscratch:10**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load CCP4**
[user@cn3144 ~]$ zcat /pdb/pdb/es/pdb3esz.ent.gz > 3esz.pdb
[user@cn3144 ~]$ zcat /pdb/pdb/f6/pdb3f6r.ent.gz > 3f6r.pdb
[user@cn3144 ~]$ superpose 3esz.pdb -s "A" 3f6r.pdb -s "A" -o 3esz_ON_3f6r.pdb

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```

CCP4 also includes an interactive graphical interface. This
application requires an [X-Windows connection](/docs/connect.html). First source the correct setup file (see above),
then type 'ccp4i' at the prompt.


![CCP4i window](/images/ccp4i.gif)




