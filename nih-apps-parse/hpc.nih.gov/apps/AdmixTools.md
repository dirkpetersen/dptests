

document.querySelector('title').textContent = 'AdmixTools on Biowulf';
AdmixTools on Biowulf


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



ADMIXTOOLS ([Patterson et al. 2012](https://genetics.med.harvard.edu/reich/Reich_Lab/Software_files/2012_Patterson_AncientAdmixture_Genetics.pdf)) 
is a software package that supports formal tests of whether admixture occurred, and makes it possible to infer admixture proportions and dates. 



The package contains 6 programs: 
* convertf: See [README.CONVERTF](/docs/AdmixTools/README.CONVERTF ) for documentation of programs for converting 
file formats.
* qp3Pop: See [README.3PopTest](/docs/AdmixTools/README.3PopTest) for details of running f\_3 test. This test can 
be used as a format test of admixture with 3 populations.
* qpBound: See [README.3PopTest](/docs/AdmixTools/README.3PopTest) for details of running qpBound. This test can 
be used for estimating bounds on the admixture proportions, given 3 
populations (2 reference and one target).
* qpDstat: See [README.Dstatistics](/docs/AdmixTools/README.Dstatistics) for details of running D-statistics. This 
is a formal test of admixture with 4 populations.
* qpF4Ratio: See [README.F4RatioTest](/docs/AdmixTools/README.F4RatioTest) for details of running F4 ratio 
estimation. This program computes the admixture proportion by taking the 
ratio of two f4 tests.
* rolloff: See [README.ROLLOFF](/docs/AdmixTools/README.ROLLOFF)/ [README.ROLLOFF\_OUTPUT](/docs/AdmixTools/README.ROLLOFF_OUTPUT) 
for details for running rolloff. This program can be used for dating admixture events.


Example files are in /usr/local/apps/AdmixTools/examples/
### References:


* Paper


Documentation
* [AdmixTools website](https://genetics.med.harvard.edu/reich/Reich_Lab/Software.html) at Harvard.
* [AdmixTools on Github](https://github.com/DReichLab/AdmixTools).


Important Notes
* Module Name: AdmixTools (see [the modules page](/apps/modules.html) for more information)
* Example files in /usr/local/apps/AdmixTools/examples/



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

[user@cn3144 ~]$ **module load AdmixTools**
[+] Loading AdmixTools 4.1 ...

[user@cn3144 examples]$ **cp -r $AMT\_DATA/\* .**
[user@cn3144 examples]$ **expfit.sh rolloff-test-parfile.par**
Running rexpfit.r 
Input file:  ASW_CEU_YRI
Output file:  expfit_ASW_CEU_YRI
Expfit logfile: expfit_ASW_CEU_YRI.flog
Output plot: expfit_ASW_CEU_YRI.pdf
 expfit.sh rolloff-test-parfile.par
Jackknife logfile: expfit_ASW_CEU_YRI.log
  
Jackknife summary: ASW_CEU_YRI.jin
Jackknife mean:          5.468
Jackknife std. err:      0.253
[user@cn3145 test]$  expfit.sh rolloff-test-parfile.par
Running rexpfit.r 
Input file:  ASW_CEU_YRI
Output file:  expfit_ASW_CEU_YRI
Expfit logfile: expfit_ASW_CEU_YRI.flog
Output plot: expfit_ASW_CEU_YRI.pdf
Jackknife logfile: expfit_ASW_CEU_YRI.log

Jackknife summary: ASW_CEU_YRI.jin
Jackknife mean:          5.468
Jackknife std. err:      0.253

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. AdmixTools.sh). For example:



```

#!/bin/bash
# this file is called admix.sh

module load AdmixTools

convertf -p par.EIGENSTRAT.PED
```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch  [--mem=#] AdmixTools.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. AdmixTools.swarm). For example:



```

rolloff -p rolloff1.par > rolloff1.log 2>&1
rolloff -p rolloff2.par > rolloff2.log 2>&1
rolloff -p rolloff3.par > rolloff3.log 2>&1

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f AdmixTools.swarm [-g #] --module AdmixTools
```

where


|  |  |  |  |
| --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | --module AdmixTools Loads the AdmixTools module for each subjob in the swarm 
 | |
 | |








