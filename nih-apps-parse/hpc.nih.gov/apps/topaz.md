

document.querySelector('title').textContent = 'topaz on Biowulf';
topaz on Biowulf


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


 Topaz is application for particle detection in cryo-electron microscopy. Topaz uses convolutional neural networks trained from positive and unlabeled examples. Topaz can also do denoising of micrographs and tomograms.. 


### References:


* Bepler, T., Morin, A., Rapp, M., Brasch, J., Shapiro, L., Noble, A.J., Berger, B.
*[Positive-unlabeled convolutional neural networks for particle picking in cryo-electron micrographs](https://pubmed.ncbi.nlm.nih.gov/31591578/)*. Nat Methods 16, 1153–1160 (2019).


Documentation
* topaz Main Site: [GitHub](https://github.com/tbepler/topaz)


Important Notes
* Module Name: topaz (see [the modules page](/apps/modules.html) for more information)
* You will need to request GPU resources to run Topaz's Denoise function (see example below)
* Test data can be found in `${TOPAZ_TEST_DATA}`



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --gres=gpu:k80:1,lscratch:50 --cpus-per-task=8 --mem=4g**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn4224 are ready for job

[user@cn4224 ~]$ **cd /lscratch/$SLURM\_JOB\_ID**
[user@cn4224 ~]$ **module load topaz/0.2.4**
[user@cn4224 ~]$ **cp ${TOPAZ\_TEST\_DATA}/topaz-tutorial-data.tar.gz .**
[user@cn4224 ~]$ **tar -xzf topaz-tutorial-data.tar.gz**

```

Preprocessing step:



```

[user@cn4224 ~]$ **mkdir -p data/EMPIAR-10025/processed**
[user@cn4224 ~]$ **mkdir -p data/EMPIAR-10025/processed/micrographs**
[user@cn4224 ~]$ **topaz preprocess -d 0 -v -s 8 -o \ 
 data/EMPIAR-10025/processed/micrographs/ \
 data/EMPIAR-10025/rawdata/micrographs/\*.mrc**
# processed: 14sep05c_c_00003gr_00014sq_00004hl_00004es_c
# processed: 14sep05c_c_00003gr_00014sq_00005hl_00003es_c
# processed: 14sep05c_c_00003gr_00014sq_00007hl_00004es_c
# processed: 14sep05c_c_00003gr_00014sq_00011hl_00003es_c
# processed: 14sep05c_c_00003gr_00015sq_00015hl_00002es_c
# processed: 14sep05c_c_00003gr_00018sq_00008hl_00003es_c
# processed: 14sep05c_c_00003gr_00018sq_00010hl_00005es_c
# processed: 14sep05c_c_00003gr_00020sq_00011hl_00002es_c
# processed: 14sep05c_c_00003gr_00020sq_00011hl_00004es_c
# processed: 14sep05c_c_00004gr_00031sq_00002hl_00002es_c
# processed: 14sep05c_c_00004gr_00031sq_00005hl_00002es_c
# processed: 14sep05c_c_00004gr_00031sq_00010hl_00002es_c
# processed: 14sep05c_c_00004gr_00032sq_00007hl_00003es_c
# processed: 14sep05c_c_00004gr_00032sq_00010hl_00003es_c
# processed: 14sep05c_c_00004gr_00032sq_00029hl_00005es_c
# processed: 14sep05c_c_00004gr_00032sq_00031hl_00002es_c
# processed: 14sep05c_c_00004gr_00032sq_00033hl_00005es_c
# processed: 14sep05c_c_00004gr_00032sq_00037hl_00002es_c
# processed: 14sep05c_c_00004gr_00032sq_00037hl_00003es_c
# processed: 14sep05c_c_00004gr_00032sq_00040hl_00002es_c
# processed: 14sep05c_c_00004gr_00032sq_00040hl_00004es_c
# processed: 14sep05c_c_00004gr_00032sq_00041hl_00005es_c
# processed: 14sep05c_c_00007gr_00013sq_00004hl_00003es_c
# processed: 14sep05c_c_00007gr_00013sq_00005hl_00002es_c
# processed: 14sep05c_c_00007gr_00013sq_00006hl_00002es_c
# processed: 14sep05c_c_00007gr_00013sq_00008hl_00003es_c
# processed: 14sep05c_c_00007gr_00013sq_00008hl_00004es_c
# processed: 14sep05c_c_00007gr_00013sq_00009hl_00002es_c
# processed: 14sep05c_c_00007gr_00013sq_00009hl_00004es_c
# processed: 14sep05c_c_00007gr_00013sq_00014hl_00004es_c

[user@cn4224 ~] **topaz convert -s 8 -o \
 data/EMPIAR-10025/processed/particles.txt \
 data/EMPIAR-10025/rawdata/particles.txt**

```

Training step:



```

[user@cn4224 ~]$ **mkdir -p saved\_models** 
[user@cn4224 ~]$ **mkdir -p saved\_models/EMPIAR-10025**
[user@cn4224 ~]$ **topaz train -n 400 \
 --num-workers=8 \
 --train-images data/EMPIAR-10025/processed/micrographs/ \
 --train-targets data/EMPIAR-10025/processed/particles.txt \
 --save-prefix=saved\_models/EMPIAR-10025/model \
 -o saved\_models/EMPIAR-10025/model\_training.txt**
# Loading model: resnet8
# Model parameters: units=32, dropout=0.0, bn=on
# Loading pretrained model: resnet8_u32
# Receptive field: 71
# Using device=0 with cuda=True
# Loaded 30 training micrographs with 1500 labeled particles
# source	split	p_observed	num_positive_regions	total_regions
# 0	train	0.00163	43500	26669790
# Specified expected number of particle per micrograph = 400.0
# With radius = 3
# Setting pi = 0.0130484716977524
# minibatch_size=256, epoch_size=1000, num_epochs=10
# Done!

```

Extraction step:



```

[user@cn4224 ~]$ **mkdir -p data/EMPIAR-10025/topaz**
[user@cn4224 ~]$ **topaz extract -r 14 -x 8 -m \
 saved\_models/EMPIAR-10025/model\_epoch10.sav \
 -o data/EMPIAR-10025/topaz/predicted\_particles\_all\_upsampled.txt \
 data/EMPIAR-10025/processed/micrographs/\*.mrc**
[user@cn4224 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. topaz.sh) similar to the following.



```

#! /bin/bash

set -e

module load topaz

topaz preprocess -d 0 -v -s 8 -o \ 
  data/EMPIAR-10025/processed/micrographs/ \
  data/EMPIAR-10025/rawdata/micrographs/*.mrc

topaz convert -s 8 -o \
  data/EMPIAR-10025/processed/particles.txt \
  data/EMPIAR-10025/rawdata/particles.txt

topaz train -n 400 \
  --num-workers=8 \
  --train-images data/EMPIAR-10025/processed/micrographs/ \
  --train-targets data/EMPIAR-10025/processed/particles.txt \
  --save-prefix=saved_models/EMPIAR-10025/model \
  -o saved_models/EMPIAR-10025/model_training.txt

topaz extract -r 14 -x 8 -m \
  saved_models/EMPIAR-10025/model_epoch10.sav \
  -o data/EMPIAR-10025/topaz/predicted_particles_all_upsampled.txt \
  data/EMPIAR-10025/processed/micrographs/*.mrc

```

Submit these jobs using the Slurm [sbatch](/docs/userguide.html) command.


Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile for the first step of the pipeline (e.g. topaz.swarm). For example:



```

topaz preprocess -d 0 -v -s 8 -o \ 
  dataset1/processed/micrographs/ \
  dataset1/rawdata/micrographs/*.mrc
topaz preprocess -d 0 -v -s 8 -o \ 
  dataset2/processed/micrographs/ \
  dataset2/rawdata/micrographs/*.mrc
topaz preprocess -d 0 -v -s 8 -o \ 
  dataset3/processed/micrographs/ \
  dataset3/rawdata/micrographs/*.mrc

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f topaz.swarm [-g #] --module topaz
```

where


|  |  |  |  |
| --- | --- | --- | --- |
| -g #  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | --module topaz  Loads the topaz module for each subjob in the swarm
 | |
 | |








