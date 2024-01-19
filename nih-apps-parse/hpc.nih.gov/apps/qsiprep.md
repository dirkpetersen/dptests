

document.querySelector('title').textContent = 'QSIPREP on Biowulf';
QSIPREP on Biowulf


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



QSIPREP configures pipelines for processing diffusion-weighted MRI (dMRI) data.



### Web site


* [Home page](https://qsiprep.readthedocs.io)
* [Github page](https://github.com/PennBBL/qsiprep)


Documentation
* [QSIPREP Documentation](https://qsiprep.readthedocs.io)


Important Notes
* qsiprep on Biowulf is a singularity container built directly from docker://pennbbl/qsiprep. However, users do not need to execute singularity directly or bind directories manually because we have provided a wrapper shell script that takes care of that. 
* Module Name: qsiprep (see [the modules page](/apps/modules.html) for more information)
* qsiprep is a scientific pipeline with the potential to overload Biowulf's central filesystem. To avoid filesystem issues we recommend the following:
	1. Limit I/O with central Biowulf directories (e.g., /data, /home) by making use of local disk (/lscratch). You can make use of lscratch to store temporary working directories by directly assigning the temporary work directory of qsiprep (--work-dir flag) to /lscratch/$SLURM\_JOB\_ID (remember to allocate enough space in /lscratch).
	* Limit memory usage. Once you determine how much memory your jobs will need (by benchmarking, but see point 5 below), It is still good idea to limit the memory used by qsiprep by using the option --mem\_mb. Something important to remember regarding memory allocation is to try to allocate some extra memory for your jobs as described in the resources described in point 6 below. 
	* Limit multi-threading. Biowulf's qsiprep has been installed as a container and some of the applications in it will hyper-thread unless you explicitly limit the number of threads. You can limit the number of threads that qsiprep is allowed to use across all processes by using the --nthreads flag.
	* Opt out of uploading qsiprep metrics to qsiprep website. You can disable the submission with the --notrack flag.
	* Make use of the flag --stop-on-first-crash flag to force a stop if issues occur.
	* Profile/benchmark qsiprep jobs: We recommend making sure a given qsiprep job can scale before launching a large number of jobs. You can do this by profiling/benchmarking small jobs (e.g., a swarm of 3 qsiprep commands), then monitoring the jobs by using either the [user dashboard](https://hpc.nih.gov/dashboard/) or the commands jobhist, sjobs, squeue, and jobload (see [biowulf utilities](https://hpc.nih.gov/docs/biowulf_tools.html)). There are resources prepared by the HPC staff that go over how to monitor your jobs for benchmarking and profiling ([video](https://youtu.be/fLMJ8-t5bm4), and [slides](https://hpc.nih.gov/training/handouts/Effective_batch_system.pdf)). Once you have profiled a swarm with a few jobs and you have refined the memory, CPU (and walltime) requirements, it is probably safe to (gradually) increase the number of qsiprep jobs. For many analyses pipelines one has no way of knowing in advance how much memory or CPUs will be actually required in an HPC environment. This is why it is very important to profile/benchmark.



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

[user@cn3144 ~]$ **module load qsiprep**

[user@cn3144 ~]$ **qsiprep -h**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. qsiprep.sh). For example:



```

#!/bin/bash
# sbatch --gres=lscratch:100 --mem=32g --cpus-per-task=48 --time=72:00:00 qsiprep.sh

set -e

function fail {
    echo "FAIL: $@" >&2
    exit 1  # signal failure
}

module load qsiprep
qsiprep /data/$USER/ds000114-download /data/$USER/qsiprep.out.ds001 \
participant --participant_label sub-02 -w /lscratch/$SLURM_JOB_ID \
--notrack --nthreads $SLURM_CPUS_PER_TASK --mem_mb $SLURM_MEM_PER_NODE \
--stop-on-first-crash --output-resolution 1.2

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--gres=lscratch:#] [--cpus-per-task=#] [--mem=#] qsiprep.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. qsiprep.swarm). For example:



```

export TMPDIR=/lscratch/$SLURM_JOB_ID; \
mkdir -p $TMPDIR/out; \
mkdir -p $TMPDIR/wrk; \
qsiprep /data/$USER/ds000114-download $TMPDIR/out \
participant --participant_label sub-01 -w $TMPDIR/wrk \
--notrack --nthreads $SLURM_CPUS_PER_TASK \
--mem_mb $SLURM_MEM_PER_NODE
--stop-on-first-crash --output-resolution 1.2; \
mv $TMPDIR/out /data/$USER/QSIPREP.out.s001
export TMPDIR=/lscratch/$SLURM_JOB_ID; \
mkdir -p $TMPDIR/out; \
mkdir -p $TMPDIR/wrk; \
qsiprep /data/$USER/ds000114-download $TMPDIR/out \
participant --participant_label sub-02 -w $TMPDIR/wrk \
--notrack --nthreads $SLURM_CPUS_PER_TASK \
--mem_mb $SLURM_MEM_PER_NODE
--stop-on-first-crash --output-resolution 1.2; \
mv $TMPDIR/out /data/$USER/QSIPREP.out.s002

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f qsiprep.swarm [--gres=lscratch:#] [-g #] [-t #] --module qsiprep
```

where


|  |  |  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- | --- | --- |
| -gres=lscratch:*#*  Number of Gigabytes of local disk space allocated per process (1 line in the swarm command file)
 | -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module qsiprep Loads the qsiprep module for each subjob in the swarm 
 | |
 | |
 | |
 | |








