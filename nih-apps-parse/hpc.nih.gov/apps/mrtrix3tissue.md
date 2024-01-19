

document.querySelector('title').textContent = 'mrtrix3tissue on Biowulf';
mrtrix3tissue on Biowulf


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



MRtrix3Tissue is an application for 3-Tissue CSD modelling and analysis (based on mrtrix).



Documentation
* [MRtrix3tissue user documentation](https://3tissue.github.io/)


Important Notes
* Module Name: mrtrix3tissue (see [the modules page](/apps/modules.html) for more information)
* Note that the GUI options have not been installed. 
* Multithreaded
* **mrtrix3tissue configuration**: Users can set up a file /home/$USER/.mrtrix.conf containing mrtrix configuration options with default options that they desire.
The list of possible settings for the mrtrix config file is [detailed here](http://mrtrix.readthedocs.io/en/latest/reference/config_file_options.html).
A sample config file is available in /usr/local/apps/mrtrix/mrtrix.conf.



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

[user@cn3144 ~]$ **module load mrtrix3tissue**
[+] Loading mrtrix3tissue  5.2.9  on cn3144

[user@cn3144 ~]$ **ss3t\_csd\_beta1 -help**
Version 3Tissue_v5.2.9           ss3t_csd_beta1
using MRtrix3 3Tissue_v5.2.9

     ss3t_csd_beta1: external MRtrix3 project

SYNOPSIS

     SS3T-CSD: beta 1 implementation

USAGE

     ss3t_csd_beta1 [ options ] in_dMRI_data in_SFWM_resp out_WM_FOD in_GM_resp
     out_GM in_CSF_resp out_CSF

        in_dMRI_data Input dMRI dataset

        in_SFWM_resp Input single-fibre WM response function text file

        out_WM_FOD   Output WM FOD image

        in_GM_resp   Input GM response function text file

        out_GM       Output GM image

        in_CSF_resp  Input CSF response function text file

        out_CSF      Output CSF image

DESCRIPTION

     This is an implementation of SS3T-CSD for beta testing and distribution.
     Use with caution and check all results carefully.

     For more information on how to use SS3T-CSD, please visit
     https://3Tissue.github.io/doc/ss3t-csd.html

OPTIONS
[...]
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. mrtrix3tissue.sh). For example:



```

#!/bin/bash
set -e
module load mrtrix3tissue

cd /data/${USER}

dwi2response dhollander dwi-prep.mif response_wm.txt response_gm.txt response_csf.txt
ss3t_csd_beta1 dwi-prep.mif response_wm.txt wmfod.mif response_gm.txt gm.mif response_csf.txt csf.mif -mask mask.mif

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] mrtrix3tissue.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. mrtrix3tissue.swarm). For example:



```

cd /data/${USER}/dataset01;\ 
   ss3t_csd_beta1 dwi-prep.mif response_wm.txt wmfod.mif response_gm.txt gm.mif response_csf.txt csf.mif -mask mask.mif
cd /data/${USER}/dataset02;\
   ss3t_csd_beta1 dwi-prep.mif response_wm.txt wmfod.mif response_gm.txt gm.mif response_csf.txt csf.mif -mask mask.mif
[...]

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f mrtrix3tissue.swarm [-g #] [-t #] --module mrtrix3tissue
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module mrtrix3tissue Loads the mrtrix3tissue module for each subjob in the swarm 
 | |
 | |
 | |








