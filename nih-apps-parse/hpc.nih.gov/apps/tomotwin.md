

document.querySelector('title').textContent = "TomoTwin on Biowulf";
TomoTwin on Biowulf


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



TomoTwin is an application the enables particle picking for Cryo-ET using deep metric learning based procedures.





> 
> TomoTwin comes pre-trained on so far 120 different proteins. By embedding tomograms in an information-rich, 
> high-dimensional space which separates macromolecules according to their 3-dimensional structure, TomoTwin allows
> users to identify proteins in tomograms de novo without manually creating training data or retraining the network
> each time a new protein is to be located. That means, you can simply run it for your specific sample without a
> much additional effort.
> 





### Reference:


* Rice, G., Wagner, T., Stabrin, M. et al.
 [**TomoTwin: generalized 3D localization of macromolecules in cryo-electron tomograms with structural data mining**](https://doi.org/10.1038/s41592-023-01878-z)
*Nat Methods (2023)*


Documentation
* [TomoTwin Site](https://tomotwin-cryoet.readthedocs.io/)
* [TomoTwin Github](https://github.com/MPI-Dortmund/tomotwin-cryoet)


Important Notes
This application uses the napari application to visualize the tomograms. napari requires a [graphical connection using NX](/docs/connect.html#nx)


* Module Name: tomotwin (see [the modules page](/apps/modules.html) for more information)
* GPU enabled
* Environment variables set 
	+ TOMOTWIN\_TEST\_DATA
	+ TOMOTWIN\_MODEL* Example files in $TOMOTWIN\_TEST\_DATA
* Model file in $TOMOTWIN\_MODEL



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive --mem=16G --gres=lscratch:50,gpu:v100x:2**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load tomotwin**

[user@cn3144 ~]$ **cd /lscratch/$SLURM\_JOB\_ID**

[user@cn3144 46116226]$ **cp -r $TOMOTWIN\_TEST\_DATA/\* .**

[user@cn3144 46116226]$ **CUDA\_VISIBLE\_DEVICES=0,1 tomotwin\_embed.py tomogram \
 -m tomotwin\_model\_p120\_052022\_loss.pth \
 -v tomo/tomo.mrc \
 -o out/embed/tomo/ \
 -b 400**
Latest version of TomoTwin is installed :-)
reading tomotwin_model_p120_052022_loss.pth
Model config:
{'identifier': 'SiameseNet', 'network_config': ...}
... UserWarning: This DataLoader will create 12 worker processes in total. Our suggested max number of worker in current system is 4 ...
Embeddings have shape: (5083356, 35)
Wrote embeddings to disk to out/embed/tomo/tomo_embeddings.temb
Done.

[user@cn3144 46116226]$ **tomotwin\_tools.py extractref \
 --tomo tomo/tomo.mrc \
 --coords ref.coords \
 --out out/extracted\_ref/**
1it [00:00, 97.91it/s]
wrote subvolume reference to out/extracted_ref/

[user@cn3144 46116226]$ **CUDA\_VISIBLE\_DEVICES=0,1 tomotwin\_embed.py subvolumes \
 -m tomotwin\_model\_p120\_052022\_loss.pth \
 -v out/extracted\_ref/reference\_0.mrc \
 -o out/embed/ref**
Latest version of TomoTwin is installed :-)
reading tomotwin_model_p120_052022_loss.pth
Model config:
{'identifier': 'SiameseNet', 'network_config': ...}
UserWarning: This DataLoader will create 12 worker processes in total. Our suggested max number of worker in current system is 4 ...
Done. Wrote results to out/embed/ref/embeddings.temb

[user@cn3144 46116226]$ **tomotwin\_map.py distance \
 -r out/embed/ref/embeddings.temb \
 -v out/embed/tomo/tomo\_embeddings.temb \
 -o out/map/**
Latest version of TomoTwin is installed :-)
Read embeddings
Map references: 100%|████████████████████████████████████| 1/1 [00:04<00:00,  4.55s/it]
Prepare output...
Wrote output to out/map/map.tmap

[user@cn3144 46116226]$ **tomotwin\_locate.py findmax -m out/map/map.tmap -o out/locate/**
Latest version of TomoTwin is installed :-)
start locate  reference_0.mrc
effective global min: 0.5
Locate class 0: 100%|███████████████████████████████| 31071/31071 [00:06<00:00, 4604.62it/s]
Call get_avg_pos
done 0
Located reference_0.mrc 847
Non-maximum-supression: 100%|███████████████████████| 847/847 [00:00<00:00, 3489.54it/s]
Particles of class reference_0.mrc: 844 (before NMS: 847)

[user@cn3144 46116226]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. tomotwin.sh). For example:



```

#!/bin/bash
set -e
cd /lscratch/$SLURM_JOB_ID
module load tomotwin
cp -r $TOMOTWIN_TEST_DATA/* .

CUDA_VISIBLE_DEVICES=0,1 tomotwin_embed.py tomogram \
    -m tomotwin_model_p120_052022_loss.pth \
    -v tomo/tomo.mrc \
    -o out/embed/tomo/ \
    -b 400
tomotwin_tools.py extractref \
    --tomo tomo/tomo.mrc \
    --coords ref.coords \
    --out out/extracted_ref/
CUDA_VISIBLE_DEVICES=0,1 tomotwin_embed.py subvolumes \
    -m tomotwin_model_p120_052022_loss.pth \
    -v out/extracted_ref/reference_0.mrc \
    -o out/embed/ref
tomotwin_map.py distance \
    -r out/embed/ref/embeddings.temb \
    -v out/embed/tomo/tomo_embeddings.temb \
    -o out/map/
tomotwin_locate.py findmax -m out/map/map.tmap -o out/locate/

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] [--gres=lscratch:#,gpu:type:2] tomotwin.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. tomotwin.swarm). For example:



```

CUDA_VISIBLE_DEVICES=0,1 tomotwin_embed.py tomogram -m model.pth -v tomo1.mrc -o out1 -b 400
CUDA_VISIBLE_DEVICES=0,1 tomotwin_embed.py tomogram -m model.pth -v tomo2.mrc -o out2 -b 400
CUDA_VISIBLE_DEVICES=0,1 tomotwin_embed.py tomogram -m model.pth -v tomo3.mrc -o out3 -b 400
CUDA_VISIBLE_DEVICES=0,1 tomotwin_embed.py tomogram -m model.pth -v tomo4.mrc -o out4 -b 400

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f tomotwin.swarm [-g #] [--gres=gpu:type:2] --module tomotwin
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | --gres=gpu:type:*2* 2 GPUs required for each process (1 line in the swarm command file). Replace type with GPU types available like v100, v100x, p100, a100, etc.
 | --module tomotwin Loads the tomotwin module for each subjob in the swarm 
 | |
 | |
 | |








