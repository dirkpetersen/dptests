

document.querySelector('title').textContent = 'ModelAngelo on Biowulf';
ModelAngelo on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
 |



ModelAngelo is an automatic atomic model building program for cryo-EM maps.



### References:


* Jamali K., Kimanius D., Scheres SHW.
 [**A Graph Neural Network Approach to Automated Model Building in Cryo-EM Maps.**](https://openreview.net/forum?id=65XDF_nwI61)
*International Conference on Learning Representations (2023).*


Documentation
* [ModelAngelo Main Site](https://github.com/3dem/model-angelo)


Important Notes
* Module Name: model-angelo (see [the modules page](/apps/modules.html) for more information)
 * Multithreaded
 * GPU dependent


Some features of ModelAngelo require the hhblits command of [hhsuite](https://hpc.nih.gov/apps/hhsuite.html).


Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive --gres=gpu:a100:1 --mem=20g -c 8**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load model-angelo**

[user@cn3144 ~]$ **model\_angelo build -v map.mrc -f sequence.fasta -o output**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. model-angelo.sh). For example:



```

#!/bin/bash
set -e
module load model-angelo
model_angelo build_no_seq -v map.mrc -o output
hhblits -i output/hmm_profiles/A.hhm -d PATH_TO_DB -o A.hhr -oa3m A.a3m -M first

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --gres=gpu:a100:1 --partition=gpu --mem=20g -c 8 model-angelo.sh
```









