

document.querySelector('title').textContent = 'ANTs on Biowulf';
ANTs on Biowulf


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



Advanced Normalization Tools (ANTs) extracts information from complex datasets that include imaging. ANTs extracts information 
from complex datasets that include imaging (Word Cloud). Paired with ANTsR (answer), ANTs is useful for managing, interpreting 
and visualizing multidimensional data. ANTs depends on the Insight ToolKit (ITK), a widely used medical image processing library 
to which ANTs developers contribute.

ANTs development is led by Brian Avants and supported by other researchers and developers at PICSL and other institutions.



Documentation
* [ANTs website](http://picsl.upenn.edu/software/ants/)


Important Notes
* Module Name: ANTs (see [the modules page](/apps/modules.html) for more information)
* Some of the ANTs scripts and executables are multi-threaded. Several of them use the variable ITK\_GLOBAL\_DEFAULT\_NUMBER\_OF\_THREADS to se the number of threads. 
* Example files in $ANTS\_TESTDATA (this variable is set when you load the module)



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=5g**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **cd /data/$USER/myruns**

[user@cn3144 ~]$ **module load ANTs**

[user@cn3144 ~]$ **cp -pr $ANTS\_TESTDATA .**

[user@cn3144 ~]$ **cd TESTDATA/stnava-asymmetry-f8ecc74** 

[user@cn3144 ~]$ **./asymmetry.sh -d 2 -a 0 -f data/symm\_t.nii.gz -m data/asymm\_s.nii.gz -o XXXX**
inputs: data/symm_t.nii.gz data/asymm_s.nii.gz XXXX 2
 CenterOfMass [133.032, 163.98]
Using double precision for computations.
Input scalar image: data/asymm_s.nii.gz
Reference image: data/asymm_s.nii.gz
[...]
 1DIAGNOSTIC,    30, 2.079582101904e-04, -4.725553011784e-04, 2.0881e+01, 5.1356e-01,
  Elapsed time (stage 0): 2.0949e+01


Total elapsed time: 2.0950e+01

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. myANTs.sh). For example:



```

#!/bin/bash
#
# this file is called myANTs.sh
#
module load ANTs
cd /data/$USER/mydir
cp -pr $ANTS_TESTDATA .
cd TESTDATA/stnava-asymmetry-f8ecc74
./asymmetry.sh -d 2 -a 0 -f data/symm_t.nii.gz -m data/asymm_s.nii.gz -o XXXX

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch myANTs.sh
```

The command above will allocate 2 CPUs and 4 GB of memory to the job, which is sufficient for this test. 


Some of the ANTs scripts and executables are multi-threaded. Several of them use the variable ITK\_GLOBAL\_DEFAULT\_NUMBER\_OF\_THREADS to set the number of threads. 
For these ANTs programs, you can add the following line in your batch script:

```

export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$SLURM_CPUS_PER_TASK

```

and submit with, for example:

```

sbatch --cpus-per-task=8 myANTs.sh

```


Other ANTs programs allow for explicit setting of the number of threads. For example, the following ANTs scripts have a '-n' option to set the number of threads:

```

antsRegistrationSpaceTime.sh
antsRegistrationSyNQuick.sh
antsRegistrationSyN.sh

```

For these scripts, you can simply add -n $SLURM\_CPUS\_PER\_TASK to the command within your batch script to set the number of threads equal to allocated CPUs. Submit with, as above, the --cpus-per-task=# sbatch parameter. 

Your job may require more than 4 GB of memory, so you may need to specify the memory when submitting. e.g.

```

$ sbatch --mem=5g myjob.bat

```

would allocate 5 GB of memory for your job. 

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. ANTs.swarm). For example:



```

cd /path/to/mydir; ./asymmetry.sh -d 2 -a 0 -f data/sym.nii.1.gz -m data/asym.nii.1.gz -o out1
cd /path/to/mydir; ./asymmetry.sh -d 2 -a 0 -f data/sym.nii.2.gz -m data/asym.nii.2.gz -o out2
cd /path/to/mydir; ./asymmetry.sh -d 2 -a 0 -f data/sym.nii.3.gz -m data/asym.nii.3.gz -o out3
[...]   

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f ANTs.swarm [-g #] [-t #] --module TEMPLATE
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g #  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t #  Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module TEMPLATE  Loads the TEMPLATE module for each subjob in the swarm 
 | |
 | |
 | |
















