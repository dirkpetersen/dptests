

document.querySelector('title').textContent = 'PETprep\_HMC: Positron Emission Tomography data preprocessing for Head Motion Correction';
**PETprep\_HMC: Positron Emission Tomography data preprocessing for Head Motion Correction** 


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



Positron Emission Tomography (PET) is a state-of-the-art neuroimaging
tool for quantification of the in vivo spatial distribution of specific
molecules in the brain. It is affected by various kinds of patient movement during a scan.
The PETprep\_HMC application allows for correction of the PET results in the presence of head motions.



### References:


* Martin Nørgaard, Melanie Ganz, Claus Svarer, Vibe G. Frokjaer, Douglas N. Greve,
Stephen C. Strother, Gitte M. Knudsen   

*Optimization of preprocessing strategies in Positron Emission Tomography
(PET) neuroimaging: A [11C]DASB PET study*
 [NeuroImage](https://www.sciencedirect.com/science/article/pii/S1053811919304471) , 199 (2019) 466–479.


Documentation
* [PETprep\_HMC Github page](https://github.com/mnoergaard/petprep_hmc)


Important Notes
* Module Name: RFdissusion (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set
	+ **PETPREP\_HMC\_HOME**  installation directory
	+ **PETPREP\_HMC\_BIN**       executable directory
	+ **PETPREP\_HMC\_SRC**       source code directory
	+ **PETPREP\_HMC\_DATA**  sample data directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive** 
[user@cn3335 ~]$ **module load petprep\_hmc** 
[+] Loading singularity  3.10.5  on cn4174
[+] Loading freesurfer  7.3.2  on cn4174
[+] Loading freetype 2.12.1 on cn4174
[+] Loading FSL 6.0.4  ...
[+] Loading petprep_hmc  0.06
[user@cn3335 ~]$
[user@cn3335 ~]$ **run.py -h**
usage: run.py [-h] --bids_dir BIDS_DIR [--output_dir OUTPUT_DIR]
              [--analysis_level {participant,group}]
              [--participant_label PARTICIPANT_LABEL [PARTICIPANT_LABEL ...]]
              [--mc_start_time MC_START_TIME] [--mc_fwhm MC_FWHM] [--mc_thresh MC_THRESH]
              [--n_procs N_PROCS] [--no_resample] [--skip_bids_validator] [-v]

BIDS App for PETPrep head motion correction workflow

optional arguments:
  -h, --help            show this help message and exit
  --bids_dir BIDS_DIR   The directory with the input dataset formatted according to the BIDS
                        standard.
  --output_dir OUTPUT_DIR
                        The directory where the output files should be stored. If you are running
                        group level analysis this folder should be prepopulated with the results
                        of theparticipant level analysis.
  --analysis_level {participant,group}
                        Level of the analysis that will be performed. Multiple participant level
                        analyses can be run independently (in parallel) using the same output_dir.
  --participant_label PARTICIPANT_LABEL [PARTICIPANT_LABEL ...]
                        The label(s) of the participant(s) that should be analyzed. The label
                        corresponds to sub- from the BIDS spec (so it does not
 include "sub-"). If this parameter is not provided all subjects should be
 analyzed. Multiple participants can be specified with a space separated
 list.
 --mc\_start\_time MC\_START\_TIME
 Start time for when to perform motion correction (subsequent frame will be
 chosen) in seconds
 --mc\_fwhm MC\_FWHM FWHM for smoothing of frames prior to estimating motion
 --mc\_thresh MC\_THRESH
 Threshold below the following percentage (0-100) of framewise ROBUST RANGE
 prior to estimating motion correction
 --n\_procs N\_PROCS Number of processors to use when running the workflow
 --no\_resample Whether or not to resample the motion corrected PET data to lowest x/y/z
 dim in original data
 --skip\_bids\_validator
 Whether or not to perform BIDS dataset validation
 -v, --version show program's version number and exit
[user@cn3335 ~]$ 


```

etc.   
   

End the interactive session:

```

[user@cn3335 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```





