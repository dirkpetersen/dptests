

document.querySelector('title').textContent = 'Gaussian on Biowulf';
Gaussian on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[GaussView](#gaussview) 
[Batch job](#sbatch) 
[Swarm of jobs](#swarm) 
[Interpreting Gaussian Errors](#errors) 
 |



Gaussian provides state-of-the-art capabilities for electronic structure modeling, and is a connected system of programs for performing semiempirical and ab initio molecular orbital (MO) calculations.



### References:


* [Official Gaussian Literature Citations](http://gaussian.com/citation/)


Documentation
* [Gaussian Main Site](http://gaussian.com)


Important Notes
* Module Name: Gaussian (see [the modules page](/apps/modules.html) for more information)
* Multi-threaded, with distributed capability using [Linda](http://gaussian.com/lindaprod/)
* Environment variables set 
	+ $g09root/$g16root
	+ $GAUSS\_SCRDIR
	+ $GAUSS\_MDEF
	+ $GAUSS\_WDEF
	+ $GAUSS\_PDEF
	+ $GAUSS\_LFLAGS
* Example files in $g09root/g09/tests/com or $g16root/g16/tests/com


Interactive use of Gaussian must be run on an interactive node with local scratch space allocated. Please see
<https://hpc.nih.gov/docs/userguide.html#int>
and
<https://hpc.nih.gov/docs/userguide.html#local>
for more information.


Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --gres=lscratch:100**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load Gaussian/G16-C01**
[user@cn3144 ~]$ **g16 < g16.com > g16.log**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```

GaussView
**GaussView** is available with the command **gview** (after loading the Gaussian module).



```
gview
```

You may see some kind of error:



```
libGL error: No matching fbConfigs or visuals found
libGL error: failed to load driver: swrast
[xcb] Unknown sequence number while processing queue
[xcb] Most likely this is a multi-threaded client and XInitThreads has not been called
[xcb] Aborting, sorry about that.
gview.exe: xcb_io.c:259: poll_for_event: Assertion `!xcb_xlib_threads_sequence_lost' failed.
Aborted
```

**NOTE:** Depending on your desktop, *gview* may require the use of NX persistent connection to display graphics. Please see <https://hpc.nih.gov/docs/connect.html> for information about NX.


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
**NOTE:** Gaussian uses local scratch disk for temporary files. Because [local disk space is not available by default, it must be deliberately allocated](http://hpc.nih.gov/docs/userguide.html#local). Otherwise you may see an error regarding local scratch space. The environment variable **$GAUSS\_SCRDIR** is set to **/lscratch/$SLURM\_JOBID** by default for biowulf jobs.


Create a batch script for submission to the cluster, for example **gaussian.batch**:



```
#---------------------  file gaussian.batch  -----------------------
#!/bin/bash
#SBATCH -J gaussian
module load Gaussian/G16-C01
g16 < g16.com > g16.log

```

The following environment variables will be set when the g16 command is run:


* **$GAUSS\_WDEF** : Run with Linda (only for multinode jobs, substitutes for %LindaWorkers), set when **--nodes** is greater than 1
* **$GAUSS\_PDEF** : Number of processors per node (substitutes for %NProcShared), set by **--cpus-per-task**
* **$GAUSS\_MDEF** : Memory per node (substitutes for %Mem), set by **--mem**


Because these environment variables are set, there is no need to include the Link 0 commands %LindaWorkers, 
%NProcShared, or %Mem. The number of processors, nodes, and amount of memory is determined by the amounts
allocated in the slurm job. However, if these commands are present in the g16 input file, they will take precedence over
the environment variables.



To submit an 8 cpu job on a single node, utilizing 4GB of memory, and 100 GB of scratch space, you would use the following command:



```
sbatch --cpus-per-task=8 --threads-per-core=1 --mem=4g --gres=lscratch:100 gaussian.batch
```

The **--threads-per-core=1** option causes slurm to allocate a single cpu per core, essentially disabling
hyperthreading, and will give better runtimes than otherwise.


Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile, for example **gaussian.swarm**:



```
#---------------------  file gaussian.swarm  -----------------------
g16 < test1.com > test1.log
g16 < test2.com > test2.log
g16 < test3.com > test3.log
g16 < test4.com > test4.log
```

Submit to swarm, allocating the required resources as necessary. The following example allocates 
4 cores (**-t 4**), 2 GB of memory (**-g 2**), and 20 GB of local disk space (**--gres=lscratch:50**)
per gaussian process.



```
swarm -f gaussian.swarm -t 4 -g 2 --noht --gres=lscratch:50 --module Gaussian
```


Interpreting Gaussian Errors
Gaussian errors are not always straightforward to interpret. Something as
simple as a "file not found" can seem baffling and cryptic. Here is a
collection of errors and their translations:




| Gaussian Error | Translation to English |
| --- | --- |
| 
```
Error termination in NtrErr:
ntran open failure returned to fopen.
Segmentation fault
```
 | Can't open a file. |
| 
```
Out-of-memory error in routine UFChkP (IEnd= 12292175 MxCore=6291456)
Use %Mem=12MW to provide the minimum amount of memory required to complete this step.
Error termination via Lnk1e at Thu Feb 2 13:05:32 2006.
```
 | Default memory (6 MW, set in $GAUSS\_MEMDEF) is too small for unfchk. |
| 
```
galloc: could not allocate memory.: Resource temporarily unavailable
```
 | Not enough memory. |
| 
```
Out-of-memory error in routine...
```
 | Not enough memory. |
| 
```
End of file in GetChg.
Error termination via Lnk1e ...
```
 | Not enough memory. |
| 
```
IMax=3 JMax=2 DiffMx= 0.00D+00
Unable to allocate space to process matrices in G2DrvN:
NAtomX= 58 NBasis= 762 NBas6D= 762 MDV1= 6291106 MinMem= 105955841.
```
 | Gaussian has 6 MW free memory (MDV1) but requires at least 106 MW
(MinMem). |
| 
```
Estimate disk for full transformation -677255533 words. Semi-Direct
transformation. Bad length for file.
```
 | MaxDisk has been set too low. |
| 
```
Error termination in NtrErr:
NtrErr Called from FileIO.
```
 | The calculation has exceeded the maximum limit of maxcyc. |
| 
```
Erroneous read. Read 0 instead of 6258688.  fd = 4 g_read
```
 | Disk quota or disk size exceeded. Could also be disk failure or NFS
timeout. |
| 
```
Erroneous write. Write 8192 instead of 12288.
fd = 4
orig len = 12288 left = 12288 g_write
```
 | Disk quota or disk size exceeded. Could also be disk failure or NFS
timeout. |
| 
```
PGFIO/stdio: Permission denied
PGFIO-F-/OPEN/unit=11/error code returned by host stdio - 13.
File name = /lscratch/2394920/Gau-#####.inp
In source file ml0.f, at line number 177
```
 | The user does not have write permission for $GAUSS\_SCRDIR. |
| 
```
QPERR — A SYNTAX ERROR WAS DETECTED IN THE INPUT LINE.
```
 | A syntax error was detected in the input. There is some guidance in the output file below this line. |
| 
```
Error during math dispatch processing...
Error: Fastmath dispatch table is corrupt
Error during math dispatch processing...
Error: Fastmath dispatch table is corrupt
```
 | g16 prior to C01 will not run on newer AMD using the AVX2 instruction set of Zen+ architecture. |




