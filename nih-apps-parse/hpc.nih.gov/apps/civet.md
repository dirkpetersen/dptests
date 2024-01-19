

document.querySelector('title').textContent = 'Civet on Biowulf';
Civet on Biowulf


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


 Civet is a pipeline for analysis of brain imaging data. Civet is used to extract and analyze cortical surfaces from MR images and can be used to perform a number of volumetric and corticometric computations. 


### References:


* C. Lepage, L. Lewis, S. Jeun, P. Bermudez, N. Khalili-Mahani, M. Omidyegaheh, A. Zijdenbos, R.D. Vincent, R. Adalat, A.C. Evans.
*[Human MR Evaluation of Cortical Thickness Using CIVET v2.1](https://archive.aievolution.com/2017/hbm1701/index.cfm?do=abs.viewAbs&abs=3292)*. Organization for Human Brain Mapping (2017)


Documentation
* Civet Main Site: [CIVET Project webpage](https://www.bic.mni.mcgill.ca/ServicesSoftware/CIVET)


Important Notes
* Module Name: civet (see [the modules page](/apps/modules.html) for more information)



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
salloc.exe: Nodes cn4224 are ready for job

[user@cn4224 ~]$ **module load civet**
[+] Loading civet  2.1.1  on cn4224 
[+] Loading singularity  3.7.1  on cn4224
[user@cn4224 ~]$ **CIVET\_Processing\_Pipeline -help**

CIVET_Processing_Pipeline, version 2.1.1, released December, 2018.

    Takes any number of multi or single spectral input MINC volumes and
    extracts the cortical surfaces from them utilizing the PMP pipeline
    system. It then calculates cortical thickness at each vertex of the 
    produced cortical surfaces (non-linearly registered) using the t-link 
    metric (in both Talairach and native spaces). It can also produce ANIMAL 
    segmentations, symmetry analyses, regional thickness, surface areas 
    and volumes for brain lobes.

    On-line documentation is available at:
        http://www.bic.mni.mcgill.ca/ServicesSoftware/CIVET.

[...]

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. civet.sh) similar to the following.



```

#! /bin/bash

set -e

module load civet

CIVET_Processing_Pipeline -sourcedir /data/user/civet-data-sets \
                          -targetdir /data/user/civet-output-dir \
                          -prefix my_data -N3-distance 0 -area-fwhm 20 \
                          -thickness tlink 30 -VBM -resample-surfaces \
                          -combine-surfaces -id-file scans -run

```

Submit these jobs using the Slurm [sbatch](/docs/userguide.html) command.


Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile for the first step of the pipeline (e.g. civet.swarm). For example:



```

CIVET_Processing_Pipeline -sourcedir /data/user/civet-data-set-01 \
                          -targetdir /data/user/civet-output-01 \
                          -prefix my_data -N3-distance 0 -area-fwhm 20 \
                          -thickness tlink 30 -VBM -resample-surfaces \
                          -combine-surfaces -id-file scans -run
CIVET_Processing_Pipeline -sourcedir /data/user/civet-data-02 \
                          -targetdir /data/user/civet-output-02 \
                          -prefix my_data -N3-distance 0 -area-fwhm 20 \
                          -thickness tlink 30 -VBM -resample-surfaces \
                          -combine-surfaces -id-file scans -run

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f civet.swarm [-g #] --module civet
```

where


|  |  |  |  |
| --- | --- | --- | --- |
| -g #  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | --module civet  Loads the civet module for each subjob in the swarm
 | |
 | |








