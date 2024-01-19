

document.querySelector('title').textContent = 'FSL on Biowulf';
FSL on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
[Swarm of jobs](#swarm) 
[FSL parallelization](#par)
[FSL on GPUs](#gpu)
 |



[![fsl logo](/images/fsl-logo-big.jpg)](http://www.fmrib.ox.ac.uk/fsl/)
FSL is a comprehensive library of image analysis and statistical tools for
FMRI, MRI and DTI brain imaging data. FSL is written mainly by members of the
Analysis Group, FMRIB, Oxford, UK. [FSL
website](http://www.fmrib.ox.ac.uk/fsl/).



FSL can analyze the following:  

Functional MRI: FEAT, MELODIC, FABBER, BASIL, VERBENA  

Structural MRI: BET, FAST, FIRST, FLIRT & FNIRT, FSLVBM, SIENA & SIENAX, fsl\_anat  

Diffusion MRI: FDT, TBSS, EDDY, TOPUP  

GLM / Stats: GLM general advice, Randomise, Cluster, FDR, Dual Regression, Mm, FLOBS  

Other: FSLView, Fslutils, Atlases, Atlasquery, SUSAN, FUGUE, MCFLIRT, Miscvis, POSSUM, BayCEST  

[Detailed Overview of FSL tools](http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslOverview)



Documentation
* [FSL website Site](http://www.fmrib.ox.ac.uk/fsl/) at Oxford University


Important Notes
* Module Name: fsl (see [the modules page](/apps/modules.html) for more information)
* All parallelization in FSL v4.0 and up is done via the fsl\_sub command that is built into several tools. fsl\_sub has been modified to use the 
Biowulf swarm utility. (Thanks to Adam Thomas and Joe Naegele NIMH)
The following programs in FSL can use parallelization: FEAT, MELODIC, TBSS,
BEDPOSTX, FSLVBM, POSSUM. See [the FSL website](http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation#Cluster_aware_tools) for
more information.

* When using FEAT you would need allocate space in lscratch and assign the TMPDIR environment variable to that allocation as follows:

```

[user@biowulf]$ export TMPDIR=${SLURM_JOB_ID}

```
* Some defaults in FSL can be overridden by setting environment variables:
	+ FSL\_MEM: Default is 4 GB, but you can increase it by setting this environment variable in a batch script or interactive session: 
	
	```
	
	[user@biowulf]$ export FSL_MEM=20
	
	```
	+ FSL\_QUEUE: By default FSL submits to the norm partition, but you can submit to another partition by setting this variable in a batch script or interactive session. Example:
	
	```
	
	[user@biowulf]$ module load fsl
	[user@biowulf]$ export FSL_QUEUE=norm
	[user@biowulf]$ bedpostx myjob
	
	```
	+ NOBATCH: Some FEAT jobs do not parallelize correctly. You can also choose to run such jobs in serial mode by setting NOBATCH=true in a batch script or interactive session. Example:
	
	```
	
	[user@biowulf]$ export NOBATCH=true
	[user@biowulf]$ bedpostx myjob
	
	```



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

[user@cn3144 ~]$ **module load fsl**

[user@cn3144 ~]$ **fsl**
![](/images/FSL_6.0.4_screenshot.png)
You should now see the FSL GUI appear on your desktop as above. Once you are finished using the GUI, please exit your interactive session by typing 'exit'.


[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$



```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. fsl.sh). For example:



```

#!/bin/bash
set -e
module load fsl/6.0.4
mcflirt -in /data/user/fmri1 -out mcf1 -mats -plots -refvol 90 -rmsrel -rmsabs; betfunc mcf1 bet1

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch  [--mem=#] fsl.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. fsl.swarm). For example:



```

mcflirt -in /data/user/fmri1 -out mcf1 -mats -plots -refvol 90 -rmsrel -rmsabs; betfunc mcf1 bet1
mcflirt -in /data/user/fmri2 -out mcf2 -mats -plots -refvol 90 -rmsrel -rmsabs; betfunc mcf2 bet2
mcflirt -in /data/user/fmri3 -out mcf3 -mats -plots -refvol 90 -rmsrel -rmsabs; betfunc mcf3 bet3
[...]

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f fsl.swarm [-g #] [-t #] --module fsl/6.0.4
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module fsl/6.0.4 Loads the fsl module for each subjob in the swarm 
 | |
 | |
 | |


FSL Parallelization
All parallelization in FSL v4.0 and up is done via the fsl\_sub command that is built into several tools. fsl\_sub has been modified to use the Biowulf swarm utility.


The following programs in FSL can use parallelization: FEAT, MELODIC, TBSS,
BEDPOSTX, FSLVBM, POSSUM. See [the FSL website](http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation#Cluster_aware_tools) for
more information.


Sample session running bedpostx in parallel: 

```

[user@biowulf]$ **module load fsl/6.0.4**

[user@biowulf]$ **bedpostx sampledataset**
subjectdir is /data/user/bedpost/sampledataset
Making bedpostx directory structure
Queuing preprocessing stages
Input args=-T 60 -m as -N bpx_preproc -l /data/user/bedpost/sampledataset.bedpostX/logs /usr/local/apps/fsl/5.0/fsl/bin/bedpostx_preproc.sh /data/user/bedpost/sampledataset 0

Queuing parallel processing stage

----- Bedpostx Monitor -----
Input args=-j 12050 -l /data/user/bedpost/sampledataset.bedpostX/logs -M user@mail.nih.gov -N bedpostx -t /data/user/bedpost/sampledataset.bedpostX/commands.txt

Queuing post processing stage
Input args=-j 12051 -T 60 -m as -N bpx_postproc -l /data/user/bedpost/sampledataset.bedpostX/logs /usr/local/apps/fsl/5.0/fsl/bin/bedpostx_postproc.sh /data/user/bedpost/sampledataset

[user@biowulf]$ 

```

The job is submitted in 3 parts: pre-processing, bedpost, and post-processing. Each part is dependent on the completion of
the previous part, and will only run after the previous part is completed. 

If you use 'sjobs' or some variant of 'squeue' to monitor your jobs, you will see at first: 

```

[user@biowulf ~]$ **sjobs**
User    JobId       JobName    Part  St  Reason      Runtime   Wallt  CPUs  Memory    Dependency          Nodelist
==================================================================================================================
user 7264_0      bpx_preproc  norm  R    ---         0:04   2:00:00      2  4GB/node                    cn1824
user 7265_[0-71] bedpostx     norm  PD   ---         0:05   2:00:00      2  4GB/node  afterany:7264_* 
user 7266_[0]    bpx_postproc norm  PD  Dependency   0:00   2:00:00      1  4GB/node  afterany:7265_*

```


ie. 3 jobs with dependencies for the 2nd and 3rd jobs.

Once the pre-processing (job 7264 in this example) is over, the main bedpost jobs (job 7265) will run, and you will see something like this:

```

[user@biowulf ~]$ **sjobs**
JOBID            NAME      TIME        ST      CPUS  MIN_ME   NODE   DEPENDENCY       NODELIST(REASON)
7264_0         bedpost     0:07        R       1     4G       1                          p1718
7264_1         bedpost     0:07        R       1     4G       1                          p1719
7264_2         bedpost     0:07        R       2     4G       1                          p999
7264_3         bedpost     0:07        R       2     4G       1                          p999
...etc...
7264_[33-71]   bedpost      0:00       PD      1     4G       1                         (QOSMaxCpusPerUserLimit)
7264_[0]       bpx_postproc 0:00       PD      1     4G       1      afterany:9517_*    (Dependency)

```


and once those are completed, the post-processing step (job 9518) will run.

#### FSL parallel jobs and memory


By default, all parallel FSL jobs submitted through fsl\_sub will be submitted requesting 4 GB memory. In some cases, this may not be enough memory for the job. In those cases, the user can
set an environment variable, FSL\_MEM, which will set the memory required for the jobs. For example, if a bedpost run as above failed due to lack of memory, you could run as follows:

```

[user@biowulf]$ **module load fsl/6.0.4**

[user@biowulf]$ **export FSL\_MEM=16**

[user@biowulf]$ **bedpostx sampledataset**
[...]

[user@biowulf]$ **sjobs**
JOBID            NAME      TIME        ST      CPUS  MIN_ME   NODE   DEPENDENCY       NODELIST(REASON)
7264_0         bedpost     0:07        R       1     16G       1                          p1718
7264_1         bedpost     0:07        R       1     16G       1                          p1719
7264_2         bedpost     0:07        R       2     16G       1                          p999
7264_3         bedpost     0:07        R       2     16G       1                          p999
...etc...
7264_[33-71]   bedpost      0:00       PD      1     16G       1                         (QOSMaxCpusPerUserLimit)
7264_[0]       bpx_postproc 0:00       PD      1     16G       1      afterany:9517_*    (Dependency)

```

In the example above, the value of FSL\_MEM was set to 16 (GB), and therefore the jobs were submitted requesting 16 GB of memory.

FSL on GPUs

Some FSL programs can use GPUs. However, they may require different versions of CUDA libraries, so there are separate modules set up for each of the GPU-enabled FSL applications. 



|  |  |  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- | --- | --- |
| **Program  **Module
| bedpostx  fsl/6.0.x/bedpostx\_gpu (note: FSL\_GPU environment variable also needs to be set)
| probtrax  fsl/6.0.x/probtrax\_gpu
| eddy  fsl/6.0.x/eddy\_cuda
 | |
 | |
 | |** |** |


See examples below for more information. 

#### bedpostx\_gpu


As of FSL 6.0.0, bedpostx\_gpu is GPU-enabled. To run bedpostx\_gpu, you need to set an additional environment variable FSL\_GPU specifying
which kind of GPU you want to run on. For example: 

```

% export FSL_GPU=k80	# submit to the K80 GPUs
or 
% export FSL_GPU=p100	# submit to the p100 GPUs
or 
% export FSL_GPU=v100	# submit to the v100 GPUs

```

The available GPU types on Biowulf can be seen by typing 'freen | grep gpu'. 

Sample session below. In 

```

[user@biowulf]$ **module load fsl/6.0.4/bedpost\_gpu**
[+] Loading fsl  CUDA/9.1  ... 
[+] Loading FSL 6.0.4  ... 

[user@biowulf]$ **export FSL\_GPU=k80**

[user@biowulf]$ **bedpostx\_gpu fdt\_subj1**
---------------------------------------------
------------ BedpostX GPU Version -----------
---------------------------------------------
subjectdir is /data/user/bedpost/fdt_subj1
Making bedpostx directory structure
Copying files to bedpost directory
Pre-processing stage
Queuing parallel processing stage

----- Bedpostx Monitor -----
Queuing post processing stage
1 parts processed out of 4
2 parts processed out of 4
3 parts processed out of 4
4 parts processed out of 4

[user@biowulf]$ **squeue -u $USER**
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
      22908533_[0]     quick bedpostx   user PD       0:00      1 (None)
      22908536_[0]     quick bedpostx   user PD       0:00      1 (Dependency)
    22908534_[0-3]       gpu bedpostx   user PD       0:00      1 (Dependency)

[user@biowulf]$ All parts processed    

```

As you see above, the pre-processing and post-processing stages run on the quick partition. Only the actual bedpostx processing is done on the GPU partition. This is transparent to the user.


#### probtrackx2\_gpu


FSL's probabilistic tracking script with crossing fibres is also GPU enabled (but not cluster-aware). probtrackx2\_gpu for CUDA 10.0 is available on Biowulf. In the example below we load the module fsl/6.0.4/probtrax\_gpu on an sinteractive session on a GPU node but for longer runs it is advisable to use sbatch and select a GPU node(s):

```

[user@biowulf]$ **sinteractive --gres=gpu:k80x:1**
salloc.exe: Pending job allocation 12345678
salloc.exe: job 12345678 queued and waiting for resources
salloc.exe: job 12345678 has been allocated resources
salloc.exe: Granted job allocation 12345678
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn1234 are ready for job

[user@cn1234]$ **module load fsl/6.0.4/probtrax\_gpu**
[+] Loading fsl CUDA/9.1 ... 
[+] Loading FSL 6.0.4  ... 

[user@cn1234]$ **cd sampledataset.bedpostX**

[user@cn1234]$ **probtrackx2\_gpu --samples=merged -m nodif\_brain\_mask.nii.gz -x nodif\_brain.nii.gz -o fdt\_paths**
PROBTRACKX2 VERSION GPU
Log directory is: logdir
Running in seedmask mode
Number of Seeds: 222441
Time Loading Data: 8 seconds
...................Allocated GPU 0...................
Free memory at the beginning: 6301548544 ---- Total memory: 6379143168
Free memory after copying masks: 5891424256 ---- Total memory: 6379143168
Running 390938 streamlines in parallel using 2 STREAMS
Total number of streamlines: 1112205000
Free memory before running iterations: 1031929856 ---- Total memory: 6379143168
Iteration 1 out of 5690
...
Iteration 5688 out of 5690
Iteration 5689 out of 5690
Iteration 5690 out of 5690

Time Spent Tracking:: 6 seconds

save results

TOTAL TIME: 14 seconds


```


#### Eddy


The eddy program can run on GPUs, via the eddy\_cuda executable. You will need to load the appropriate module so that eddy\_cuda can
find the appropriate CUDA libraries. In FSL 6.0.5, there
are two different compiled versions that you can use on Biowulf: eddy\_cuda9.1 and eddy\_cuda10.2. Sample batch script: 


```

#/bin/bash 

cd mydir
module load fsl/6.0.5/eddy_cuda
eddy_cuda10.2 [...]

```


You will also need to submit to a GPU node, of course. e.g.

```

sbatch --partition=gpu --gres=gpu:k80:1 myjob.sh

```





























