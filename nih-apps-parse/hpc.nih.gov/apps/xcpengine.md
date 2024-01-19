

document.querySelector('title').textContent = 'xcpEngine on Biowulf';
xcpEngine on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
 |



xcpEngine is a pipeline application that performs denoising of fMRI datasets and computes functional connectivity. This pipeline can take as input a dataset that has been preprocessed with [fmriprep](https://hpc.nih.gov/apps/fmriprep.html)



### Web site


* [Home page](https://xcpengine.readthedocs.io)
* [Github page](https://github.com/PennBBL/xcpEngine)


Documentation
* [xcpEngine Documentation](https://xcpengine.readthedocs.io)


Important Notes
* Module Name: xcpengine (see [the modules page](/apps/modules.html) for more information)
* xcpEngine is a scientific pipeline with the potential to overload Biowulf's central filesystem. To avoid filesystem issues we recommend the following:
	1. Limit I/O with central Biowulf directories (e.g., /data, /home) by making use of local disk (/lscratch). You can make use of lscratch to store temporary working directories by directly assigning the temporary work directory of xcpEngine (-i option) to /lscratch/$SLURM\_JOB\_ID (remember to allocate enough space in /lscratch). You should also use re-direct the output of xcpEngine to /lscratch/$SLURM\_JOB\_ID (-o option), then copy the output back to your data directory at the end of your job (see [batch example](https://hpc.nih.gov/apps/xcpengine.html#sbatch) below). Additionally, you shuold copy your input dataset and cohort and design files to /lscratch/$SLURM\_JOB\_ID to keep I/O local as much as possible (see [batch example](https://hpc.nih.gov/apps/xcpengine.html#sbatch) below).
	* Profile/benchmark xcpEngine jobs: We recommend making sure a given xcpEngine job can scale before launching a large number of jobs. You can do this by profiling/benchmarking small jobs (e.g., a swarm of 3 xcpEngine commands), then monitoring the jobs by using either the [user dashboard](https://hpcnihapps.cit.nih.gov/auth/dashboard/) or the commands jobhist, sjobs, squeue, and jobload (see [biowulf utilities](https://hpc.nih.gov/docs/biowulf_tools.html)). There are resources prepared by the HPC staff that go over how to monitor your jobs for benchmarking and profiling ([video](https://youtu.be/fLMJ8-t5bm4), and [slides](https://hpc.nih.gov/training/handouts/Effective_batch_system.pdf)). Once you have profiled a swarm with a few jobs and you have refined the memory, CPU (and walltime) requirements, it is probably safe to (gradually) increase the number of xcpEngine jobs. For many analyses pipelines one has no way of knowing in advance how much memory or CPUs will be actually required in an HPC environment. This is why it is very important to profile/benchmark.



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session to display usage (user input in **bold**):



```

[user@biowulf]$ **sinteractive**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn1234 are ready for job

[user@cn1234 ~]$ **module load xcpengine**

[user@cn1234 ~]$ **xcpEngine**
   Usage: xcpEngine -d  

 Compulsory arguments:
 -d : Primary design file for pipeline:
 The design file specifies the pipeline modules to
 be included in the current analysis and configures
 any parameters necessary to run the modules.
 -c : Cohort file for pipeline input:
 A comma-separated catalogue of the analytic sample.
 Each row corresponds to a subject, and each column
 corresponds either to an identifier or to an input.
 -o : Parent directory for pipeline output:
 A valid path on the current filesystem specifying
 the directory wherein all output from the current
 analysis will be written.

 Optional arguments:
 -m : Execution mode:
 Input can either be 's' (for serial execution on a
 single machine)[default], 'c' (for execution on a
 computing cluster) or a path to a file (for execution
 on a computing cluster, subject to the specifications
 defined in the file).
 -i : Scratch space for pipeline intermediates:
 Some systems operate more quickly when temporary
 files are written in a dedicated scratch space. This
 argument enables a scratch space for intermediates.
 -r : Root directory for inputs:
 If all paths defined in the cohort file are defined
 relative to a root directory, then this argument will
 define the root directory. Otherwise, all paths will
 be treated as absolute.
 -t : Trace:
 Integer value ( 0 - 3 ) that indicates the level
 of verbosity during module execution. Higher
 levels reduce readability but provide useful
 information for troubleshooting.
 0 [default]: Human-readable explanations of
 processing steps and error traces only.
 1: Explicitly trace module-level computations.
 Print a workflow map when execution completes.
 2: Explicitly trace module- and utility-level
 computations.
 3: All commands called by the module and all
 children are traced and explicitly replicated
 in a log file.

[user@cn1234 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. xcpengine.sh). For example:



```

#!/bin/bash
# sbatch --gres=lscratch:100 --mem=32g --cpus-per-task=48 --time=72:00:00 xcpengine.sh

set -e

function fail {
    echo "FAIL: $@" >&2
    exit 1  # signal failure
}

module load xcpengine; \
cd /lscratch/$SLURM_JOB_ID/
cp /data/user/XCPENGINE_TEST/fmriprep.tar.gz /lscratch/$SLURM_JOB_ID/.; \
tar -xzf /lscratch/$SLURM_JOB_ID/fmriprep.tar.gz -C /lscratch/$SLURM_JOB_ID; \
cp /data/user/XCPENGINE_TEST/anat-antsct.dsn /lscratch/$SLURM_JOB_ID/.; \
cp /data/user/XCPENGINE_TEST/anat_cohort.csv /lscratch/$SLURM_JOB_ID/.; \
xcpEngine -d /lscratch/$SLURM_JOB_ID/anat-antsct.dsn -c /lscratch/$SLURM_JOB_ID/anat_cohort.csv -o /lscratch/$SLURM_JOB_ID/xcp_output -t 1 -i /lscratch/$SLURM_JOB_ID -r /lscratch/$SLURM_JOB_ID; \
tar -czf xcp_output.tar.gz xcp_output; \
cp /lscratch/$SLURM_JOB_ID/xcp_output.tar.gz /data/user/XCPENGINE_TEST/. 

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--gres=lscratch:#] [--cpus-per-task=#] [--mem=#] xcpengine.sh
```







