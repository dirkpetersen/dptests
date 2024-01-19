

document.querySelector('title').textContent = 'XCP-D on Biowulf';
XCP-D on Biowulf


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



XCP-D uses the output of FMRIPREP to generate denoised BOLD images, parcellated time-series, functional connectivity matrices, and quality assesment reports.



### References:


* Adebimpe, Azeez, Bertolero, Maxwell, Mehta, Kahini, Salo, Taylor, Murtha, Kristin, Cieslak, Matthew, Meisler, Steven, Madison, Thomas, Sydnor, Valerie, Covitz, Sydney, Fair, Damien, & Satterthwaite, Theodore.
 [*XCP-D : A Robust Postprocessing Pipeline of fMRI data.*](https://zenodo.org/record/7717239) 
 Zenodo. https://doi.org/10.5281/zenodo.7717239.


  

Documentation
* [XCP-D Documentation](https://xcp-d.readthedocs.io)
* [Github page](https://github.com/PennLINC/xcp_d)


Important Notes
* Module Name: xcp\_d (see [the modules page](/apps/modules.html) for more information)
* Test data in $XCP\_D\_TEST\_DATA



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

[user@cn3144 ~]$ **module load xcp\_d**
[+] Loading xcp_d  0.5.0  on cn3144 
[+] Loading singularity  3.10.5  on cn3144 

[user@cn3144 ~]$ **xcp\_d --help**
usage: xcp_d [-h] [--version] [--participant_label PARTICIPANT_LABEL [PARTICIPANT_LABEL ...]]
[-t TASK_ID] [-m] [-s] [--nthreads NTHREADS] [--omp-nthreads OMP_NTHREADS]
[--mem_gb MEM_GB] [--use-plugin USE_PLUGIN] [-v] [--input-type {fmirprep,dcan,hpc}]
[--smoothing SMOOTHING] [--despike]
[-p {27P,36P,24P,acompcor,aroma,acompcor_gsr,aroma_gsr}]
[-c CUSTOM_CONF] [-d DUMMYTIME] [--lower-bpf LOWER_BPF] [--upper-bpf UPPER_BPF]
[--bpf-order BPF_ORDER] [--motion-filter-type {lp,notch}]
[--band-stop-min BAND_STOP_MIN] [--band-stop-max BAND_STOP_MAX]
[--motion-filter-order MOTION_FILTER_ORDER] [-r HEAD_RADIUS]
[-f FD_THRESH] [-w WORK_DIR] [--clean-workdir] [--resource-monitor] [--notrack]
fmri_dir output_dir

xcp_d postprocessing workflow of fMRI data

positional arguments:
  fmri_dir              the root folder of a preprocessed fMRI output .
  output_dir            the output path for xcp_d

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit

[...]

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. xcp\_d.sh). For example (using test data below):



```

#!/bin/bash
#SBATCH --job-name=xcp_d
#SBATCH --gres=lscratch:20
#SBATCH --mem=8g
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=4

module load xcp_d/0.5.0

tar -C /lscratch/${SLURM_JOB_ID} -xf ${XCP_D_TEST_DATA}/fmriprep.out.ds001.tar.gz

xcp_d /lscratch/${SLURM_JOB_ID}/fmriprep.out.ds001 /data/$USER/xcp_d.out\
         participant --participant_label sub-01 -w /lscratch/${SLURM_JOB_ID}  

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch xcp_d.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. xcp\_d.swarm). For example:



```

xcp_d /data/${USER}/BIDS-dataset/fmriprep.out.ds001/ /data/$USER/XCP-D/xcp_d.out.ds001 \
         participant --participant_label sub-01 -w /lscratch/${SLURM_JOB_ID} \
xcp_d /data/${USER}/BIDS-dataset/fmriprep.out.ds002/ /data/$USER/XCP-D/xcp_d.out.ds002 \
         participant --participant_label sub-02 -w /lscratch/${SLURM_JOB_ID} \
xcp_d /data/${USER}/BIDS-dataset/fmriprep.out.ds003/ /data/$USER/XCP-D/xcp_d.out.ds003 \
         participant --participant_label sub-03 -w /lscratch/${SLURM_JOB_ID} \

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f xcp_d.swarm [--gres=lscratch:#] [-g #] -t auto --module xcp_d
```

where


|  |  |  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- | --- | --- |
| -gres=lscratch:*#*  Number of Gigabytes of local disk space allocated per process (1 line in the swarm command file)
 | -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file). We set this to `auto` to allocate all CPUs in each node.
 | --module xcp\_d Loads the xcp\_d module for each subjob in the swarm 
 | |
 | |
 | |
 | |








