

document.querySelector('title').textContent = 'DeepLoc2 on Biowulf';
DeepLoc2 on Biowulf


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



 DeepLoc2 uses deep learning to predict subcellular localization of eukaryotic proteins. 
 
> 
>  DeepLoc 2.0 predicts the subcellular localization(s) of eukaryotic proteins. 
>  DeepLoc 2.0 is a multi-label predictor, which means that is able to predict one or 
>  more localizations for any given protein. It can differentiate between 10 different 
>  localizations: Nucleus, Cytoplasm, Extracellular, Mitochondrion, Cell membrane, 
>  Endoplasmic reticulum, Chloroplast, Golgi apparatus, Lysosome/Vacuole and Peroxisome. 
>  Additionally, DeepLoc 2.0 can predict the presence of the sorting signal(s) that had 
>  an influence on the prediction of the subcellular localization(s).
>  





### References:


* Vineet Thumuluri, Jose Juan Almagro Armenteros, Alexander Rosenberg Johansen, Henrik Nielsen, Ole Winther.
 [**DeepLoc 2.0: multi-label subcellular localization prediction using protein language models.**](https://doi.org/10.1093/nar/gkac278)
*Nucleic Acids Research, Web server issue 2022.*


Documentation
* [DeepLoc2 Main Site](https://services.healthtech.dtu.dk/service.php?DeepLoc-2.0)


Important Notes
* Module Name: deeploc (see [the modules page](/apps/modules.html) for more information)
* singlethreaded
* Environment variables set 
	+ DEEPLOC\_TRAIN\_DATA
	+ DEEPLOC\_TEST\_DATA* Example files in $DEEPLOC\_TEST\_DATA or /usr/local/apps/deeploc/TEST\_DATA
* The first time you run deeploc2, it will **download checkpoint data**. Due to the size of the checkpoint
 files, deeploc2 on Biowulf is configured to download this to /data/$USER/.cache/torch/hub/checkpoints.
 Users can also copy the files from $DEEPLOC\_TRAIN\_DATA/checkpoints before running deeploc2 to skip this step. 
 However, for this to work you can only copy to **/data/$USER/.cache/torch/hub/checkpoints**.
* **Memory Considerations:** We recommend at least 24G memory for your runs. Please run test jobs to benchmark your dataset.



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program on the test data. Then compare it with the test data result using diff.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive --mem=24G --gres=lscratch:5**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load deeploc**

[user@cn3144 ~]$ **mkdir -p /data/$USER/.cache/torch/hub**

[user@cn3144 ~]$ **cp -r $DEEPLOC\_TRAIN\_DATA/checkpoints /data/$USER/.cache/torch/hub/**

[user@cn3144 ~]$ **cd /lscratch/$SLURM\_JOB\_ID**

[user@cn3144 ~]$ **cp $DEEPLOC\_TEST\_DATA/test.fasta .**

[user@cn3144 ~]$ **deeploc2 -f test.fasta**

[user@cn3144 ~]$ **diff outputs/results\_20230101-000000.csv $DEEPLOC\_TEST\_DATA/results\_test.csv**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. deeploc.sh). For example:



```

#!/bin/bash
set -e
module load deeploc
cd /data/$USER
deeploc2 -f input.fasta

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--mem=#] deeploc.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. deeploc.swarm). For example:



```

deeploc2 -f 01.fasta -o results_01
deeploc2 -f 02.fasta -o results_02
deeploc2 -f 03.fasta -o results_03
deeploc2 -f 04.fasta -o results_04

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f deeploc.swarm [-g #] --module deeploc
```

where


|  |  |  |  |
| --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | --module deeploc Loads the deeploc module for each subjob in the swarm 
 | |
 | |








