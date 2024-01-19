

document.querySelector('title').textContent = 'rosettafold on Biowulf';
rosettafold on Biowulf


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



Accurate prediction of protein structures and interactions using a three-track neural network, in which information at the 1D sequence level, the 2D distance map level, and the 3D coordinate level is successively transformed and integrated. 






### References:


* Baek M et al.*Accurate prediction of protein structures and interactions using a three-track neural network* Science. 2021 Jul 15
 [PubMed](https://pubmed.ncbi.nlm.nih.gov/34282049/) | 
 [Journal](https://science.sciencemag.org/content/early/2021/07/19/science.abj8754.long)


Documentation
* RoseTTAFold Github:[Github](https://github.com/RosettaCommons/RoseTTAFold)


Important Notes
* Module Name: RoseTTAFold (see [the modules page](/apps/modules.html) for more information)
 * RoseTTAFold needs write permission to the model, so please save a local copy of network at your home directory after you load the module the first time.
 
```

	module load RoseTTAFold
	cp -r ${ROSETTAFOLD_NETWORK:-none} ~/
	
```
* To run complex modeling, you also needs to save a local copy of weights at your home directory. 
 
```

	cp -r ${ROSETTAFOLD_WEIGHTS:-none} ~/
	
```
* For PPI screening using faster 2-track version (only available for RoseTTAFold/1.1.0), you need to copy a different network to home directory.
 
```

	cp -r ${ROSETTAFOLD_NETWORK_2TRACK:-none} ~/
	
```
* Since only the last step of run\_e2e\_ver.sh and run\_pyrosetta\_ver.sh can use GPU, we strongly suggest to run the edited pipeline which was spited to part1(CPUs and memory heavy) and part2(GPU), see the examples in interactive job.




Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive --cpus-per-task=10 --mem=60G**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144]$ **module load RoseTTAFold**
[user@cn3144]$ **mkdir /data/$USER/rosettafold\_test/**
[user@cn3144]$ **cd /data/$USER/rosettafold\_test/**
[user@cn3144]$ **cp -r ${ROSETTAFOLD\_TEST\_DATA:-none}/\* .**
[user@cn3144]$ **run\_e2e\_ver\_part1.sh input.fa e2e\_out**
Running HHblits
Running PSIPRED
Running hhsearch
Running end-to-end prediction
Done with part1, please run part2 on GPU node

[user@cn3144]$ **run\_pyrosetta\_ver\_part1.sh input.fa pyrosetta\_out**
Running HHblits
Running PSIPRED
Running hhsearch
Predicting distance and orientations
Running parallel RosettaTR.py
Done with part1, please run part2 at GPU node

[user@cn3144 ]$ **exit**
salloc.exe: Relinquishing job allocation 46116226

[user@biowulf]$ **sinteractive --cpus-per-task=2 --mem=10g --gres=gpu:p100:1**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144]$ **module load RoseTTAFold**
[user@cn3144]$ **cd /data/$USER/rosettafold\_test/**
[user@cn3144]$ **run\_e2e\_ver\_part2.sh input.fa e2e\_out**
run_e2e_ver_part2.sh input.fa e2e_out
Running end-to-end prediction
Done with part2 (prediction)

[user@cn3144]$ **run\_pyrosetta\_ver\_part2.sh input.fa pyrosetta\_out**
Picking final models
Final models saved in: pyrosetta_out/model
Done with part2 (pick final models)

```

 For PPI screening using faster 2-track version:

```

[user@biowulf]$ **sinteractive --cpus-per-task=2 --mem=10g --gres=gpu:p100:1**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144]$ **module load RoseTTAFold**
[user@cn3144]$ **mkdir /data/$USER/rosettafold\_test/**
[user@cn3144]$ **cd /data/$USER/rosettafold\_test/**
[user@cn3144]$ **cp -r ${ROSETTAFOLD\_TEST\_DATA:-none}/\* .**
[user@cn3144]$ **cd complex\_2track**
[user@cn3144]$ **python ~/network\_2track/predict\_msa.py -msa input.a3m -npz complex\_2track.npz -L1 218**

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. rosettafold.sh). For example:



```


#!/bin/bash
set -e
module load RoseTTAFold
cd /data/$USER/rosettafold_test/
cp -r ${ROSETTAFOLD_TEST_DATA:-none}/* .
cd complex_modeling
python ~/network/predict_complex.py -i paired.a3m -o complex3 -Ls 218 310


```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=2 --mem=10g --partition=gpu --gres=gpu:v100x:1 rosettafold.sh
```











