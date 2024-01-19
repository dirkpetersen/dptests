

document.querySelector('title').textContent = 'Acemd on Biowulf';
Acemd on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
[Swarm of jobs](#swarm) 
[Benchmarks](/apps/acemd/benchmarks.html)
 |



HTMD is a molecular-specific programmable environment to prepare, handle, simulate, visualize and analyse molecular systems. HTMD includes the Acemd software. 



HTMD/Acemd is produced by Acellera. [[Acellera website](https://www.acellera.com)]

The Acemd license has expired, and usage of the program on Biowulf is not sufficient
to justify the cost of a license. Acemd will run only in 'Basic' mode, which limits you to a single GPU. If you need 
higher performance on multiple GPUs, please migrate your work to [using
other MD programs on Biowulf](https://hpc.nih.gov/apps/#compchem) instead of Acemd

Documentation
* [HTMD/Acemd Documentation](https://software.acellera.com/acemd/index.html)


Important Notes
* Module Name: acemd (see [the modules page](/apps/modules.html) for more information)
* Acemd is designed to run on GPUs. It will run on the the p100s, but is not yet supported on the v100 GPUs.
* There are 256 floating Acemd licenses available on Biowulf by courtesy of NIAID. You can check the status of the acemd licenses by typing 'licenses'.
* Sample files for Acemd3 are in /usr/local/apps/acemd/acemd3\_examples/



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --gres=gpu:p100:1** 
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load acemd** 
[+] Loading acemd  3.5  ...

[user@cn3144 ~]$  **cd mydir**

[user@cn3144 ~]$ **cp -r $(dirname $(which acemd3))/../share/acemd3/dhfr\_charmm .**
[user@cn3144 ~]$ **cd dhfr\_charmm** 
[user@cn3144 ~]$ **acemd3 input**
#
# ACEMD version 3.5
#
# Copyright (C) 2017-2019 Acellera (www.acellera.com)
#
# When publishing, please cite:
#   ACEMD: Accelerating Biomolecular Dynamics in the Microsecond Time Scale
#   M. J. Harvey, G. Giupponi and G. De Fabritiis,
#   J Chem. Theory. Comput. 2009 5(6), pp1632-1639
#   DOI: 10.1021/ct9000685
#
# Licence:
#   Check floating licence:
# Licence:
#   Check floating licence:
#     ACELLERA_LICENCE_SERVER -- not defined
#     ACELLERA_LICENSE_SERVER -- not defined
#   Check node-locked licence:
#     ACELLERA_LICENCE_FILE -- not defined
#     ACELLERA_LICENSE_FILE -- not defined
#     /opt/acellera/licence.dat -- DENIED (Unable to locate target file)
#     /opt/acellera/license.dat -- DENIED (Unable to locate target file)
#     /home/user/.acellera/licence.dat -- DENIED (Unable to locate target file)
#     /home/user/.acellera/license.dat -- DENIED (Unable to locate target file)
#
# ACEMD is running with a basic licence!
#
# Contact Acellera (info@acellera.com) for licensing.
#
# WARNING: ACEMD is limited to run on the GPU device 0 only!
# Read input file: input
# Parse input file
[...]
# Step       Time         Bond         Angle        Urey-Bradley Dihedral     Improper     CMAP         Non-bonded   Implicit
 External     Potential    Kinetic      Total        Temperature  Volume
#            [ps]         [kcal/mol]   [kcal/mol]   [kcal/mol]   [kcal/mol]   [kcal/mol]   [kcal/mol]   [kcal/mol]   [kcal/mol]
 [kcal/mol]   [kcal/mol]   [kcal/mol]   [kcal/mol]   [K]          [A^3]
       25000       100.00     537.7867    1281.4603     155.2703    1634.3376      95.5368     -66.3912  -74851.6306       0.0000
       0.0000  -71213.6300   14366.8914  -56846.7386      298.867    240990.21
# Speed: average  204.99 ns/day, current  204.99 ns/day
# Progress: 0.4%, remaining time: 3:15:59, ETA: Thu May  7 18:03:31 2020
[...]
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. Acemd.sh). For example:



```

#!/bin/bash
set -e
cd /data/$USER/somedir
module load acemd
acemd3 input

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --partition=gpu --gres=gpu:p100:1  jobscript
```


Benchmarks

Please see the [Acemd benchmarks page](/apps/acemd/benchmarks.html).
















