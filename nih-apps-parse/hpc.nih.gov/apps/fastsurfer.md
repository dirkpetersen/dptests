

document.querySelector('title').textContent = 'Fastsurfer on Biowulf';
Fastsurfer on Biowulf


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



Fastsurfer is a neuroimaging pipeline based on deep learning. 



References
* Henschel L, Conjeti S, Estrada S, Diers K, Fischl B, Reuter M. [FastSurfer - A fast and accurate deep learning based neuroimaging pipeline](https://pubmed.ncbi.nlm.nih.gov/32526386/), Neuroimage. 2020;219:117012


Documentation
* [Fastsurfer github page](https://github.com/Deep-MI/FastSurfer)


Important Notes
* Module Name: fastsurfer (see [the modules page](/apps/modules.html) for more information). Note that you have to load the freesurfer and python modules as well. The freesurfer initialization scripts need to be run in addition to loading the module. 

* Freesurfer uses temporary files for some calculations. By default, these will go to /scratch (i.e. the top level of the shared scratch directory) which is low-performance (i.e. reading and writing to /scratch is slowing down your jobs). Instead, you should allocate local disk (/lscratch) and set the environment variable $tmpdir=/lscratch/${SLURM\_JOBID} to speed up your jobs. See the examples below.



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) with 5 GB of local disk and run the program. Sample session below::



```

[user@biowulf]$ **sinteractive --mem=35g --cpus-per-task=64 --gres=lscratch:5**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load fastsurfer freesurfer python**
[+] Loading fastsurfer  c5e9677  on cn3144 
[+] Loading freesurfer  7.1.1  on cn3144 
[+] Loading python 3.7  ...

[user@cn3144 ~]$ **cd /lscratch/${SLURM\_JOB\_ID}**

[user@cn3144 ~]$ **source $FREESURFER\_HOME/SetUpFreeSurfer.sh**
-------- freesurfer-linux-centos7_x86_64-7.1.1-20200723-8b40551 --------
Setting up environment for FreeSurfer/FS-FAST (and FSL)
FREESURFER_HOME   /usr/local/apps/freesurfer/7.1.1
FSFAST_HOME       /usr/local/apps/freesurfer/7.1.1/fsfast
FSF_OUTPUT_FORMAT nii.gz
SUBJECTS_DIR      /usr/local/apps/freesurfer/7.1.1/subjects
MNI_DIR           /usr/local/apps/freesurfer/7.1.1/mni

[user@cn3144 ~]$ **export tmpdir=/lscratch/${SLURM\_JOBID}**

[user@cn3144 ~]$ **run\_fastsurfer.sh --t1 $SUBJECTS\_DIR/bert/mri/orig.mgz \
 --sid bert \
 --sd /lscratch/${SLURM\_JOB\_ID}/analysis \
 --parallel --threads 64**
Thu Aug  6 15:17:30 EDT 2020

/usr/local/apps/fastsurfer/c5e9677/FastSurferCNN /lscratch/46116226
python eval.py --in_name /usr/local/apps/freesurfer/7.1.1/subjects/bert/mri/orig.mgz --out_name /lscratch/46116226/analysis/bert/mri/aparc.DKTatlas+aseg.deep.mgz --order 1 --network_sagittal_path ../checkpoints/Sagittal_Weights_FastSurferCNN/ckpts/Epoch_30_training_state.pkl --network_axial_path ../checkpoints/Axial_Weights_FastSurferCNN/ckpts/Epoch_30_training_state.pkl --network_coronal_path ../checkpoints/Coronal_Weights_FastSurferCNN/ckpts/Epoch_30_training_state.pkl --batch_size 8 --simple_run
Reading volume /usr/local/apps/freesurfer/7.1.1/subjects/bert/mri/orig.mgz
Loading Axial
Successfully loaded Image from /usr/local/apps/freesurfer/7.1.1/subjects/bert/mri/orig.mgz
Loading Sagittal
Successfully loaded Image from /usr/local/apps/freesurfer/7.1.1/subjects/bert/mri/orig.mgz
Loading Coronal.
Successfully loaded Image from /usr/local/apps/freesurfer/7.1.1/subjects/bert/mri/orig.mgz

[...]

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. fastsurfer.sh). For example:



```

#!/bin/bash
set -e
module load fastsurfer freesurfer python
source $FREESURFER_HOME/SetUpFreeSurfer.sh
# set the environment variable tmpdir to local scratch for better performance
export tmpdir=/lscratch/$SLURM_JOBID

run_fastsurfer.sh --t1 $SUBJECTS_DIR/bert/mri/orig.mgz \
                  --sid bert \
                  --sd /lscratch/${SLURM_JOB_ID}/analysis \
                  --parallel --threads 64

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch  [--mem=#] --gres=lscratch:5 fastsurfer.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
For swarm jobs, it might be simplest to have the line

```

module load fastsurfer freesurfer python > /dev/null 2>&1 ; source $FREESURFER_HOME/SetUpFreeSurfer.sh

```

in a bash script file. Alternatively, you can add this line to *each line* in your swarm command file.


Create a swarmfile (e.g. fastsurfer.swarm). For example:



```

export tmpdir=/lscratch/$SLURM_JOBID;run_fastsurfer.sh --t1 $SUBJECTS_DIR/bert/mri/orig.mgz --sid bert --sd /lscratch/${SLURM_JOB_ID}/analysis
export tmpdir=/lscratch/$SLURM_JOBID;run_fastsurfer.sh --t1 $SUBJECTS_DIR/bert/mri/orig.mgz --sid bert --sd /lscratch/${SLURM_JOB_ID}/analysis 
export tmpdir=/lscratch/$SLURM_JOBID;run_fastsurfer.sh --t1 $SUBJECTS_DIR/bert/mri/orig.mgz --sid bert --sd /lscratch/${SLURM_JOB_ID}/analysis 

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f fastsurfer.swarm [-t #] --gres=lscratch:5
```

where


|  |  |  |  |
| --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | --gres=lscratch:5 allocate 5 GB of local disk for each swarm subjob
 | |
 | |










