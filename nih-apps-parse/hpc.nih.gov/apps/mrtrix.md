

document.querySelector('title').textContent = 'mrtrix on Biowulf';
mrtrix on Biowulf


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



MRtrix provides a large suite of tools for image processing, analysis and visualisation, with a focus on the analysis of white matter using diffusion-weighted MRI



Documentation
* [mrtrix user documentation](http://mrtrix.readthedocs.io/en/latest/index.html)


Important Notes
* Module Name: mrtrix (see [the modules page](/apps/modules.html) for more information)
* Note that the GUI options have not been installed. 
* Multithreaded
* Environment variables set 
	+ MRTRIX\_HOME* **mrtrix configuration**: Users can set up a file /home/$USER/.mrtrix.conf containing mrtrix configuration options with default options that they desire.
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

[user@cn3144 ~]$ **module load mrtrix**
[+] Loading eigen 3.3.4  ...
[+] Loading Zlib 1.2.11  ...
[+] Loading GSL 2.4 for GCC 4.8.5 ...
[+] Loading FFTW 3.3.7 , compiled with gcc4.8.5  and openmpi2.1.2  ...
[+] Loading mrtrix 3.0_RC2  ...

[user@cn3144 ~]$ **mrinfo /usr/local/apps/mrtrix/sample\_data/data\_slice\_0000.nii.gz**
mrinfo: [WARNING] transform matrix contains invalid entries - resetting to sane defaults
************************************************
Image:               "/usr/local/apps/mrtrix/sample_data/data_slice_0000.nii.gz"
************************************************
  Dimensions:        128 x 104 x 1 x 65
  Voxel size:        2 x 2 x 2 x 1
  Data strides:      [ 1 2 4 3 ]
  Format:            NIfTI-1.1 (GZip compressed)
  Data type:         signed 16 bit integer (little endian)
  Intensity scaling: offset = 0, multiplier = 1
  Transform:                    1           0           0        -127
                                0           1           0        -103
                                0           0           1          -0
  comments:          FSL5.0
  
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. mrtrix.sh). For example:



```

#!/bin/bash
set -e
module load mrtrix
mrstats /usr/local/apps/mrtrix/sample_data/nodif_brain_mask.nii.gz
mrhistogram /usr/local/apps/mrtrix/sample_data/nodif.nii.gz  out.hist

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] mrtrix.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. mrtrix.swarm). For example:



```

mrstats file1.nii.gz
mrstats file2.nii.gz
mrstats file3.nii.gz
[...]

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f mrtrix.swarm [-g #] [-t #] --module mrtrix
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module mrtrix Loads the mrtrix module for each subjob in the swarm 
 | |
 | |
 | |








