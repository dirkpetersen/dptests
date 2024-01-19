

document.querySelector('title').textContent = 'Solar on Biowulf';
Solar on Biowulf


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



SOLAR is a package of software to perform several kinds of statistical
genetic analysis, including linkage analysis, quantitative genetic analysis,
and covariate screening. The name SOLAR stands for "Sequential Oligogenic
Linkage Analysis Routines."


Solar was developed by researchers at the Southwest Foundation for
Biomedical Research. 
[Solar website](http://solar-eclipse-genetics.org/index.html).


Documentation
* [Solar website](http://solar-eclipse-genetics.org/index.html)
* [solar commands](http://solar-eclipse-genetics.org/support.html)
* [Solar FAQ](http://solar-eclipse-genetics.org/faq.html)


Important Notes
* Module Name: solar (see [the modules page](/apps/modules.html) for more information)
* Singlethreaded, GPU-capabale
* Example files in /usr/local/apps/solar/examples


Users looking to take advantage of new GPU capabilities should load the solar/8.5.1b or solar/9.0.1 module and allocate a GPU node as necessary



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

[user@cn3144 ~]$ **cp -r $SOLAR\_EXAMPLES .** 

[user@cn3144 ~]$  **cd examples**

[user@cn3144 ~]$ $ **solar**

SOLAR Eclipse version 9.0.0, last updated on February 04, 2022
Developed at Maryland Psychiatric Research Center,
University of Maryland School of Medicine, Baltimore.
Visit our documentation and tutorial website www.solar-eclipse-genetics.org
Our download page https://www.nitrc.org/projects/se_linux
Our github page https://github.com/brian09/solar-eclipse
For questions email: pkochunov@gmail.com
Enter help for help, exit to exit, doc to browse documentation.
The software development is supported by NIH grant RO1EB015611
from The National Institute for Biomedical Imaging and Bioengineering.
Enter cite to see how to cite this software.

solar> **load pedigree gaw10.ped**
Unloading current pedigree data ...
Loading pedigree data from the file gaw10.ped ...
solar> **load phenotypes phen**
phen: ID AGE EF Q1 Q2 Q3 Q4 Q5
solar> **trait Q4**
solar> **covariate sex age age\*sex**
solar> **polygenic**
**********************************************************************
*  Maximize sporadic model                                           *
**********************************************************************

    *** Parameters    e2   h2r
    ***       Values  1 0
    *** Upper Bounds  1.01 1
    *** Lower Bounds  0.03 -0.0001

                     Building index for phenotypes file...

IARRAY allocation: 64646
RARRAY allocation: 32930
                              Pedigree:  gaw10.ped
                               Phenotypes:  phen

                     Merging pedigree and phenotype data...

                          No pedigrees were skipped.

                              Pedigrees Included
                              --------- --------
           Pedigrees  People   Females   Males    MZTwins Probands
               23      1497      756      741        0        0
[...etc...]
******************************************************************************
*                          Summary of Results                                *
******************************************************************************

        Pedigree:    gaw10.ped
        Phenotypes:  phen
        Trait:       Q4                    Individuals:  1000

                         H2r is 0.5518983  p = 2.8983051e-29  (Significant)
               H2r Std. Error:  0.0575342


        Proportion of Variance Due to All Final Covariates Is
                                  0.0080729

        Output files and models are in directory Q4/
        Summary results are in Q4/polygenic.out
        Loglikelihoods and chi's are in Q4/polygenic.logs.out
        Best model is named poly and null0 (currently loaded)
        Final models are named poly, spor, nocovar

        Residual Kurtosis is -0.0786, within normal range
solar> **quit**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. solar.sh). For example:



```

#!/bin/bash
set -e
module load solar

cp /usr/local/apps/solar/examples/* .

solar << EOF
load pedigree gaw10.ped
load phenotypes phen
trait Q4
covariate sex age age*sex
polygenic
quit
EOF

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] solar.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
1. First create different directories for each solar run. Put all the
required input files under the created directories.


2. For each directory, create a script file which contains the solar
commands as below:





```

-----------/data/$USER/solar/run1/script ----------
module load solar 
cd /data/$USER/solar/run1
solar << EOF
load pedigree ped.txt
load marker mrk.txt
verbosity min
freq mle
load map map.txt
.......
....
...
return ""
quit
EOF
-------------------------------------------------

```

3. Now prepare the swarm command file with one line for each Solar run,
e.g.




```

#------- cmdfile -------------
/data/$USER/solar/run1/script
/data/$USER/solar/run2/script
/data/$USER/solar/run3/script
/data/$USER/solar/run4/script
.....
....
/data/$USER/solar/runX/script
#---- end of cmdfile ---------

```

Submit this swarm of Solar jobs with

```

swarm -f cmdfile  --g #] --module solar

```

where 


|  |  |  |  |
| --- | --- | --- | --- |
| -f cmdfile   swarm command file
| -g #  # of Gigabytes of memory required for a single Solar run
 | |
 | |










