

document.querySelector('title').textContent = 'SIDESPLITTER on Biowulf';
SIDESPLITTER on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
 |



Sidesplitter reduces over-fitting in both idealised and experimental settings, while maintaining
independence between the two sides of a split refinement. It can improve the final resolution in refinements of structures prone to severe over-fitting, such as membrane proteins in detergent micelles.



Documentation
* [SIDESPLITTER Main Site](https://github.com/StructuralBiology-ICLMedicine/SIDESPLITTER)


Important Notes
* Module Name: sidesplitter (see [the modules page](/apps/modules.html) for more information)
 * Multithreaded (**NOTE:** sidesplitter will attempt to use all cpus on a node by default unless constrained by **OMP\_NUM\_THREADS**)



Sidesplitter can also be used in tandom with [RELION](/apps/RELION) v3.1:



```

# Usage:
#     sidesplitter_wrapper.sh path/to/relion_external_reconstruct_star_file.star
#
# To use from RELION 3.1, set the RELION_EXTERNAL_RECONSTRUCT_EXECUTABLE environment variable to point to this script,
# set SIDESPLITTER to point to the sidesplitter binary (or make sure sidesplitter can be found via PATH), and run
# relion_refine with the --external_reconstruct argument. For example:
#
#     export RELION_EXTERNAL_RECONSTRUCT_EXECUTABLE=/path/to/sidesplitter_wrapper.sh
#     export SIDESPLITTER=/path/to/sidesplitter
#
# then run RELION auto-refine from the GUI and put "--external_reconstruct" in the additional arguments box. To run on
# a cluster, depending on your configuration you might need to put the environment variable definitions into your
# submission script.

# Troubleshooting
#
# If you have problems running SIDESPLITTER using this script, the first thing to check is that external reconstruction
# from RELION is working correctly. Try running a normal refinement job, using the "--external_reconstruct" argument
# but without setting the RELION_EXTERNAL_RECONSTRUCT_EXECUTABLE environment variable. If this fails, the problem is
# likely to be with your RELION installation - perhaps it is the wrong version, or different installations are
# conflicting with each other. If normal external reconstruction is successful, the problem is likely to be with the
# SIDESPLITTER installation, or a bug in this script.

# How this script works:
#
# If the target file name contains "_half", this script assumes two copies of itself will be running (for the two half
# data sets). The two copies will coordinate with each other by creating, checking for and deleting a directory called
# "sidesplitter_running". Both scripts will run relion_external_reconstruct for their given half data set. The first
# script will wait for both reconstructions to finish, then call SIDESPLITTER to process both half maps. The second
# script will wait for the first to finish running SIDESPLITTER and then exit (because if either of the scripts exits
# before the processing is finished, RELION moves on and tries to continue its own processing before the filtered
# volumes are ready).
#
# If the target file name does not contain "_half", this script assumes there is only a single copy of itself running.
# In this case it call relion_external_reconstruct, waits for the reconstruction to finish and then exits.
# This handles the final iteration when the two half sets are combined, at which point RELION calls the external
# reconstruction program just once to reconstruct the final combined volume.
#
# Note that this script is not particularly robust. If one of the commands fails, it's possible you might need to
# manually tidy up the job directory and remove the "sidesplitter_running" directory to avoid problems in the next run.

```

For use with the GUI as detailed above, before launching the GUI



```
module load RELION sidesplitter
...
export RELION_EXTERNAL_RECONSTRUCT_EXECUTABLE=sidesplitter_wrapper.sh
relion
```

Then within the Running tab for a refinement run:


![Additional arguments](sidesplitter_RELION.png)
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

[user@cn3144 ~]$ **module load sidesplitter**

[user@cn3144 ~]$ **OMP\_NUM\_THREADS=$SLURM\_CPUS\_ON\_NODE sidesplitter --v1 half\_map1.mrc --v2 half\_map2.mrc --mask mask.mrc**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. sidesplitter.sh). For example:



```

#!/bin/bash
set -e
module load sidesplitter
OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE sidesplitter --v1 half_map1.mrc --v2 half_map2.mrc --mask mask.mrc

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] sidesplitter.sh
```







