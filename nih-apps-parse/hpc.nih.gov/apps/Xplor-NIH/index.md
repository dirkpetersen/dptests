

document.querySelector('title').textContent = 'Xplor-NIH on Biowulf';
Xplor-NIH on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Batch job](#sbatch) 
 |



XPLOR-NIH is a structure determination program which builds on the X-PLOR v
3.851 program, including additional tools developed at the NIH. These tools
include functionality for the following:


* 3J couplings
* 1J couplings
* 13C shifts
* 1H shifts
* T1/T2
* dipolar couplings
* radius of gyration
* CSA
* conformational database torsion angle potentials
* database base-base positioning potentials for DNA
* interface to the NMR graphics package [VMD-XPLOR](http://vmd-xplor.cit.nih.gov/), downloadable separately.
* embedded Python and TCL interpreters.







### References:


[XPLOR](http://www.csb.yale.edu/userguides/datamanip/xplor/xplorman/htmlman.html) was developed at Yale
University, and the NIH extensions were developed by G. Marius Clore, John
Kuszewski, Charles D. Schwieters, and Nico Tjandra at the NIH.


Documentation
* [XPLOR-NIH web site](http://bit.niddk.nih.gov/xplor-nih/),
including [tools for use at
NIH](http://bit.niddk.nih.gov/xplor-nih/nih/)



- [Xplor Manual](http://bit.niddk.nih.gov/xplor-nih/xplorMan/)

Important Notes
	* Module Name: Xplor-NIH (see [the modules page](/apps/modules.html) for more information)
	* Additional environment variables set: XPLOR\_NIH\_TOPPAR = location of topology/parameter files for that version
	* Xplor-NIH requires the usual set of Xplor input files -- coordinates, xplor
	input files, topology/parameter files. The 'standard' files are in
	$XPLOR\_NIH\_TOPPAR, the envronment variable that is set when the module is loaded. 
	* To run a job on Biowulf, please use the slurmXplor helper.
	* Sometimes a large Xplor-NIH slurm job will hang at startup time. In this case, try adding the option -startup\_delay 1 to the command-line.

Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).

```

 
[user@biowulf ~]$ slurmXplor -py -ntasks  num -partition partitionName [options] script.py

```

where num is the number of processes to run and options are any
 option to the xplor command. It is also possible to specify
 arbitrary sbatch options. For details please see 



```

 

[user@biowulf ~]$ slurmXplor -help

```
