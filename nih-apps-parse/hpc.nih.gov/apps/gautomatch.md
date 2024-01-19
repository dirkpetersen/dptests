

document.querySelector('title').textContent = 'Gautomatch on Biowulf';
Gautomatch on Biowulf


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



Fully automatic acccurate, convenient and extremely fast particle picking for EM 



### References:


* http://www.mrc-lmb.cam.ac.uk/kzhang/Gautomatch


Documentation
* [Gautomatch Main Site](https://www.mrc-lmb.cam.ac.uk/kzhang/Gautomatch/)


Important Notes
* Module Name: gautomatch (see [the modules page](/apps/modules.html) for more information)
* GPU app
* Example files in /usr/local/apps/gautomatch/TEST\_DATA/



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --gres=gpu:k80:1 --mem=20g -c14**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ cp /usr/local/apps/gautomatch/TEST_DATA/ice_carbon_aggregation_low-contrast/* .

[user@cn3144 ~]$ module load gautomatch
[+] Loading gautomatch 0.56 on cn4178
[+] Loading CUDA Toolkit 8.0.44 ...

[user@cn3144 ~]$ gautomatch --apixM 1.34 --diameter 400 --T templates_lp40_3.2A.mrcs --apixT 3.2 --lave_D 100 --lave_min -0.8 --lsigma_cutoff 1.2  --cc_cutoff 0.25 test?.mrc
***************************************************************************************************
User input parameters:
       --apixM              1.340
       --diameter           400.00
       --T                  templates_lp40_3.2A.mrcs
       --apixT              3.200
       --lave_D             100.000
       --lave_min           -0.800
       --lsigma_cutoff      1.200
       --cc_cutoff          0.250
[...]
File test2.mrc finished.
>>>>>>TIME<<<<<<                         PICKING: 1.049951s
#################################################################################################################################
All 2 files finished successfully:
>>>>>>TIME<<<<<<                           TOTAL: 4.786598s
#################################################################################################################################
[...]

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. gautomatch.sh). For example:



```

#!/bin/bash
set -e
module load gautomatch
gautomatch --apixM 1.34 --diameter 400 --T templates_lp40_3.2A.mrcs --apixT 3.2 --lave_D 100 --lave_min -0.8 --lsigma_cutoff 1.2  --cc_cutoff 0.25 test?.mrc

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --partition=gpu --gres=gpu:k80:1 --cpus-per-task=14 --mem=20g gautomatch.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. gautomatch.swarm). For example:



```

gautomatch [...] test1_?.mrc
gautomatch [...] test2_?.mrc
gautomatch [...] test3_?.mrc

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f gautomatch.swarm -g 20 -t 14 --module gautomatch
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module gautomatch Loads the gautomatch module for each subjob in the swarm 
 | |
 | |
 | |








