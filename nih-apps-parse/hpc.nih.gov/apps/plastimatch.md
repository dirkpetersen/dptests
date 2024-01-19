

document.querySelector('title').textContent = 'Plastimatch on Biowulf';
Plastimatch on Biowulf


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


 Plastimatch is an application to perform registration of medical images such as X-Rays, CT, PET, and MRI.. 


### References:


* Zaffino P, Raudaschl P, Fritscher K, Sharp GC, Spadea MF.
*[Technical Note: plastimatch mabs, an open source tool for automatic image segmentation](https://pubmed.ncbi.nlm.nih.gov/27587045/)*. Med Phys. 2016 Sep;43(9):5155.


Documentation
* Plastimatch Main Site: [Plastimatch](https://plastimatch.org)


Important Notes
* Module Name: plastimatch (see [the modules page](/apps/modules.html) for more information)
* Test data can be found in `${PLASTIMATCH_TEST_DATA}`



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --cpus-per-task=2 --mem=4g**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn4224 are ready for job

[user@cn4224 ~]$ **module load plastimatch**
[user@cn4244 ~]$ **cd /data/${USER}**
[user@cn4224 ~]$ **cp ${PLASTIMATCH\_TEST\_DATA}/\* .**
[user@cn4224 ~]$ **tar xvzf registration-tutorial.tar.gz**
[user@cn4224 ~]$ **cd registration-tutorial**
[user@cn4224 ~]$ **plastimatch register parms.txt**
Loading fixed image [0]: t5.mha
Loading moving image [0]: t0.mha
Launching registration worker thread
Inside registration worker thread
Doing registration stage
[1] xf_in->m_type = 0, xf_out->m_type = 0
RESAMPLE 0 1: (3 3 2), (3 3 2)
RESAMPLE 0 1: (3 3 2), (3 3 2)
volume_calc_grad complete.
plm_warp_native is complete.
[...]
Saving image...
Trying to write image to warped.mha
Load:   0.138638
Run:    1.27459
Save:   0.392006
Total:  1.80523
Finished!


```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. plastimatch.sh) similar to the following.



```

#! /bin/bash

module load plastimatch
cd /data/${USER}
cp ${PLASTIMATCH_TESTDATA}/* .
tar xzf registration-tutorial.tar.gz
cd registration-tutorial
plastimatch register parms.txt

```

Submit these jobs using the Slurm [sbatch](/docs/userguide.html) command.


Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile for the first step of the pipeline (e.g. plastimatch.swarm). For example:



```

plastimatch register parms_01.txt
plastimatch register parms_02.txt
plastimatch register parms_03.txt

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f plastimatch.swarm [-g #] --module plastimatch
```

where


|  |  |  |  |
| --- | --- | --- | --- |
| -g #  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | --module plastimatch  Loads the plastimatch module for each subjob in the swarm
 | |
 | |








