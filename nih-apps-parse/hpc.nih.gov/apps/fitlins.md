

document.querySelector('title').textContent = 'fitlins on Biowulf';
fitlins on Biowulf


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



Fitlins is an application for estimating linear models for a given BIDS-compliant dataset.



### Web site


* [Home page](https://fitlins.readthedocs.io)
* [Github page](https://github.com/poldracklab/fitlins)


  

Documentation
* [fitlins Documentation](https://fitlins.readthedocs.io)


Important Notes
* fitlins on Biowulf is a singularity container built directly from docker://poldracklab/fitlins. However, users do not need to execute singularity directly or bind directories manually because we have provided a wrapper shell script that takes care of that. 
* Module Name: fitlins (see [the modules page](/apps/modules.html) for more information)
* fitlins is a scientific pipeline with the potential to overload Biowulf's central filesystem. To avoid filesystem issues we recommend the following:
	1. Limit I/O with central Biowulf directories (e.g., /data, /home) by making use of local disk (/lscratch). You can make use of lscratch to store temporary working directories by directly assigning the temporary work directory of fitlins (-w option) to /lscratch/$SLURM\_JOB\_ID (remember to allocate enough space in /lscratch).
	* Limit memory usage. Once you determine how much memory your jobs will need (by benchmarking, but see point 5 below), It is still good idea to limit the memory used by fitlins by using the option --mem\_gb. Something important to remember regarding memory allocation is to try to allocate some extra memory for your jobs as described in the resources described in point 6 below. 
	* Limit multi-threading. Biowulf's fitlins has been installed as a container and some of the applications in it will hyper-thread unless you explicitly limit the number of threads. You can limit the number of threads that fitlins is allowed to use across all processes by using the --n-cpus flag.
	* Profile/benchmark fitlins jobs: We recommend making sure a given fitlins job can scale before launching a large number of jobs. You can do this by profiling/benchmarking small jobs (e.g., a swarm of 3 fitlins commands), then monitoring the jobs by using either the [user dashboard](https://hpc.nih.gov/dashboard/) or the commands jobhist, sjobs, squeue, and jobload (see [biowulf utilities](https://hpc.nih.gov/docs/biowulf_tools.html)). There are resources prepared by the HPC staff that go over how to monitor your jobs for benchmarking and profiling ([video](https://youtu.be/fLMJ8-t5bm4), and [slides](https://hpc.nih.gov/training/handouts/Effective_batch_system.pdf)). Once you have profiled a swarm with a few jobs and you have refined the memory, CPU (and walltime) requirements, it is probably safe to (gradually) increase the number of fitlins jobs. For many analyses pipelines one has no way of knowing in advance how much memory or CPUs will be actually required in an HPC environment. This is why it is very important to profile/benchmark.



Interactive job
In the helix session below we fetch datasets from Poldrack's Stanford Lab using datalad then process them with fitlins, as shown in the fitlins docs linked above. If you are already have your own dataset and model you can skip to the next step.   
Sample session (user input in **bold**):



```

[user@helix ~]$ **module load datalad git**
[+] Loading datalad  0.13.0rc2  on helix 
[+] Loading singularity  3.6.4  on helix
[+] Loading git 2.29.2  ... 

[user@helix ~]$ **cd /data/$USER**

[user@helix ~]$ **datalad install -r ///labs/poldrack/ds003\_fmriprep**
install(ok): /data/user/ds003_fmriprep (dataset)                                                                        
[INFO   ] Installing Dataset(/data/user/ds003_fmriprep) to get /data/user/ds003_fmriprep recursively    
[INFO   ] Remote origin not usable by git-annex; setting annex-ignore                                                                   
[INFO   ] access to 1 dataset sibling s3-PRIVATE not auto-enabled, enable with:                                                         
| 		datalad siblings -d "/data/user/ds003_fmriprep/sourcedata" enable -s s3-PRIVATE 
[INFO   ] Submodule HEAD got detached. Resetting branch master to point to 571d8737. Original location was c0905372 
install(ok): /data/user/ds003_fmriprep/sourcedata (dataset)
action summary:
  install (ok: 2)

[user@helix ~]$ **cd ds003\_fmriprep**

[user@helix ~]$ **datalad get ds003\_fmriprep/sub-\*/func/\*\_space-MNI152NLin2009cAsym\_desc-\*.nii.gz \
 ds003\_fmriprep/sub-\*/func/\*\_desc-confounds\_\*.tsv**
get(ok): sub-05/func/sub-05_task-rhymejudgment_space-MNI152NLin2009cAsym_desc-aseg_dseg.nii.gz (file) [from origin...]                  
get(ok): sub-07/func/sub-07_task-rhymejudgment_space-MNI152NLin2009cAsym_desc-aparcaseg_dseg.nii.gz (file) [from origin...]             
get(ok): sub-05/func/sub-05_task-rhymejudgment_space-MNI152NLin2009cAsym_desc-aparcaseg_dseg.nii.gz (file) [from origin...]
get(ok): sub-03/func/sub-03_task-rhymejudgment_space-MNI152NLin2009cAsym_desc-brain_mask.nii.gz (file) [from origin...]
get(ok): sub-04/func/sub-04_task-rhymejudgment_space-MNI152NLin2009cAsym_desc-brain_mask.nii.gz (file) [from origin...]
get(ok): sub-10/func/sub-10_task-rhymejudgment_space-MNI152NLin2009cAsym_desc-smoothAROMAnonaggr_bold.nii.gz (file) [from origin...]
get(ok): sub-12/func/sub-12_task-rhymejudgment_desc-confounds_regressors.tsv (file) [from origin...]
get(ok): sub-12/func/sub-12_task-rhymejudgment_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz (file) [from origin...]
get(ok): sub-05/func/sub-05_task-rhymejudgment_desc-confounds_regressors.tsv (file) [from origin...]
get(ok): sub-04/func/sub-04_task-rhymejudgment_space-MNI152NLin2009cAsym_desc-smoothAROMAnonaggr_bold.nii.gz (file) [from origin...]
  [68 similar messages have been suppressed]
action summary:
  get (ok: 78)

[user@helix ~]$ **exit**

```

[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. In the example below we fetch datasets from Poldrack's Stanford Lab using datalad then process them with fitlins, as shown in the docs linked above.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load fitlins**

[user@cn3144 ~]$ **cd /data/$USER**

[user@cn3144 ~]$ **fitlins ds003\_fmriprep/sourcedata output/ dataset \
 --derivatives $PWD/ds003\_fmriprep \
 --model model.json \
 --smoothing 5:run \
 -w /lscratch/${SLURM\_JOB\_ID} \
 --n-cpus 2**


```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. fitlins.sh). For example:



```

#!/bin/bash
# sbatch --gres=lscratch:50 --mem=32g --cpus-per-task=48 --time=06:00:00 fitlins.sh

set -o pipefail
set -e

function fail {
    echo "FAIL: $@" >&2
    exit 1  # signal failure
}

module load fitlins
fitlins /data/$USER/ds003_fmriprep/sourcedata /data/$USER/output/ dataset \
                    --derivatives /data/$USER/ds003_fmriprep \
                    --model /data/$USER/model.json \
                    --smoothing 5:run \
                    -w /lscratch/${SLURM_JOB_ID} \
                    --n-cpus 2

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--gres=lscratch:#] [--cpus-per-task=#] [--mem=#] fitlins.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. fitlins.swarm). For example:



```

export TMPDIR=/lscratch/$SLURM_JOB_ID; \
mkdir -p $TMPDIR/wrk; \
fitlins /data/$USER/ds003_fmriprep/sourcedata /data/$USER/output3/ dataset \
                    --derivatives /data/$USER/ds003_fmriprep \
                    --model /data/$USER/model.json \
                    --smoothing 5:run \
                    -w ${TMPDIR}/wrk \
                    --n-cpus 2
export TMPDIR=/lscratch/$SLURM_JOB_ID; \
mkdir -p $TMPDIR/wrk; \
fitlins /data/$USER/ds002_fmriprep/sourcedata /data/$USER/output2/ dataset \
                    --derivatives /data/$USER/ds002_fmriprep \
                    --model /data/$USER/model.json \
                    --smoothing 5:run \
                    -w ${TMPDIR}/wrk \
                    --n-cpus 2
export TMPDIR=/lscratch/$SLURM_JOB_ID; \
mkdir -p $TMPDIR/wrk; \
fitlins /data/$USER/ds001_fmriprep/sourcedata /data/$USER/output1/ dataset \
                    --derivatives /data/$USER/ds001_fmriprep \
                    --model /data/$USER/model.json \
                    --smoothing 5:run \
                    -w ${TMPDIR}/wrk \
                    --n-cpus 2

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f fitlins.swarm [--gres=lscratch:#] [-g #] [-t #] --module fitlins
```

where


|  |  |  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- | --- | --- |
| -gres=lscratch:*#*  Number of Gigabytes of local disk space allocated per process (1 line in the swarm command file)
 | -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module fitlins Loads the fitlins module for each subjob in the swarm 
 | |
 | |
 | |
 | |








