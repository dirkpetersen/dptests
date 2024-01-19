

document.querySelector('title').textContent = 'dcm2niix on Biowulf';
dcm2niix on Biowulf


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


 dcm2niix is a designed to convert neuroimaging data from the DICOM format
to the NIfTI format. 


Documentation
* dcm2niix [Manual](https://www.nitrc.org/plugins/mwiki/index.php/dcm2nii:MainPage)
* dcm2niix [GitHub](https://github.com/rordenlab/dcm2niix)


Important Notes
* Module Name: dcm2niix (see [the modules page](/apps/modules.html) for more information)
* The developer of dcm2niix, Chris Rorden, has provided the Biowulf group with the following comments to help HPC users of dmc2nix:
	+ If you want nii.gz compressed output, note that you have two options:
		- "-z y" will write the raw data to disk and then compress with the parallel pigz (assuming pigz is found)
		 - "-z i" will use the internal zip compressor: this writes uses a single thread to write compressed data directly to disk.
	In general, "-z y" is faster for single user computers that have fast SSD drives. However, clusters often have slow disk IO and the user may only want to allocate a single thread so that other tasks can run concurrently. For this reason, it may be faster to use "-z i" in a cluster situation.



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load dcm2niix/1.0.20171017**
[user@cn3144 ~]$ **dcm2niix -h**

Chris Rorden's dcm2niiX version v1.0.20171017 (OpenJPEG build) GCC7.2.0 (64-bit Linux)
usage: dcm2niix [options] 
 Options :
 -1..-9 : gz compression level (1=fastest..9=smallest, default 6)
 -b : BIDS sidecar (y/n/o(o=only: no NIfTI), default y)
 -ba : anonymize BIDS (y/n, default y)
 -c : comment stored as NIfTI aux\_file (up to 24 characters)
 -d : diffusion volumes sorted by b-value (y/n, default n)
 -f : filename (%a=antenna (coil) number, %c=comments, %d=description, %e echo number, %f=folder name, %i ID of patient, %j seriesInstanceUID, %k studyInstanceUID, %m=manufacturer, %n=name of patient, %p=protocol, %s=series number, %t=time, %u=acquisition number, %v=vendor, %x=study ID; %z sequence name; default '%f\_%p\_%t\_%s')
 -h : show help
 -i : ignore derived, localizer and 2D images (y/n, default n)
 -m : merge 2D slices from same series regardless of study time, echo, coil, orientation, etc. (y/n, default n)
 -o : output directory (omit to save to input folder)
 -p : Philips precise float (not display) scaling (y/n, default y)
 -s : single file mode, do not convert other images in folder (y/n, default n)
 -t : text notes includes private patient details (y/n, default n)
 -v : verbose (n/y or 0/1/2 [no, yes, logorrheic], default 0)
 -x : crop (y/n, default n)
 -z : gz compress images (y/i/n/3, default n) [y=pigz, i=internal, n=no, 3=no,3D]
 Defaults file : /home/user/.dcm2nii.ini
 Examples :
 dcm2niix /Users/chris/dir
 dcm2niix -c "my comment" /Users/chris/dir
 dcm2niix -o /users/cr/outdir/ -z y ~/dicomdir
 dcm2niix -f %p\_%s -b y -ba n ~/dicomdir
 dcm2niix -f mystudy%s ~/dicomdir
 dcm2niix -o "~/dir with spaces/dir" ~/dicomdir
Example output filename: '/myFolder\_MPRAGE\_19770703150928\_1.nii'

[user@cn3144 ~]$ **cp $DCM2NIIX\_TEST\_DATA/CR-MONO1-10-chest.gz .**
[user@cn3144 ~]$ **gunzip CR-MONO1-10-chest.gz**
[user@cn3144 ~]$ **dcm2niix /CR-MONO1-10-chest**
Chris Rorden's dcm2niiX version v1.0.20171017 (OpenJPEG build) GCC7.2.0 (64-bit Linux)
Found 1 DICOM image(s)
Convert 1 DICOM as /home/user/notes/software/dcm2niix/dcm2niix\_\_0\_1a (440x440x1x1)
Warning: Check that 2D images are not mirrored.
Conversion required 0.047990 seconds (0.020000 for core code).

# or for a whole folder
[user@cn3144 ~]$ **dcm2niix -o ./my\_nifty -z y ./my\_dicom**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. dcm2niix.sh), which uses the input file 'dcm2niix.in'. For example:



```

#!/bin/bash
module load dcm2niix/1.0.20171017
dcm2niix  -o ./my_nifty -z y ./my_dicom

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] dcm2niix.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. dcm2niix.swarm). For example:



```

dcm2niix  -o ./my_nifty1 -z y ./my_dicom1
dcm2niix  -o ./my_nifty2 -z y ./my_dicom2
dcm2niix  -o ./my_nifty3 -z y ./my_dicom3

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f dcm2niix.swarm [-g #] [-t #] --module dcm2niix
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g #  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t #  Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module dcm2niix  Loads the dcm2niix module for each subjob in the swarm 
 | |
 | |
 | |








