

document.querySelector('title').textContent = 'rosettafoldna on Biowulf';
rosettafoldna on Biowulf


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



rosettafoldnaNA: rapidly produces 3D structure models with confidence estimates for protein-DNA and protein-RNA complexes, and for RNA tertiary structures.






### References:


* Baek M et al.*Accurate prediction of nucleic acid and protein-nucleic acid complexes using rosettafoldnaNA* biorxiv
 [biorxiv](https://www.biorxiv.org/content/10.1101/2022.09.09.507333v1) |


Documentation
* rosettafoldna Github:[Github](https://github.com/uw-ipd/rosettafoldna2NA)


Important Notes
* Module Name: rosettafoldna (see [the modules page](/apps/modules.html) for more information)
 * Since only the last step of rosettafoldna can use GPU, we strongly suggest to run the edited pipeline which was spited to part1(CPUs and memory heavy) and part2(GPU >V100), see the examples in interactive job.




Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive --cpus-per-task=10 --mem=70G**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144]$ **module load rosettafoldna**
[user@cn3144]$ **mkdir /data/$USER/rosettafoldna\_test/**
[user@cn3144]$ **cd /data/$USER/rosettafoldna\_test/**
[user@cn3144]$ **cp -r ${ROSETTAFOLDNA\_TEST\_DATA:-none}/\* .**
[user@cn3144]$ **tree rosettafoldna\_test/**
rosettafoldna_test/
├── protein.fa
└── RNA.fa

0 directories, 2 files

[user@cn3144]$ **run\_RF2NA\_part1.sh test\_o protein.fa R:RNA.fa**
Running HHblits
Running PSIPRED
Running hhsearch
Running rMSA (lite)
Done with part1, please run part2 on GPU node (>= V100)

[user@cn3144]$ **exit**
salloc.exe: Relinquishing job allocation 46116226

```


```

[user@biowulf]$ **sinteractive --cpus-per-task=2 --mem=10g --gres=gpu:v100:1**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job
[user@cn3144]$ **run\_RF2NA\_part2.sh test\_o protein.fa R:RNA.fa**
Running RoseTTAFold2NA to predict structures
Running on GPU
  msa[msa == "U"] = 30
           plddt    best
RECYCLE  0   0.874  -1.000
RECYCLE  1   0.892   0.874
RECYCLE  2   0.898   0.892
RECYCLE  3   0.899   0.898
RECYCLE  4   0.901   0.899
RECYCLE  5   0.902   0.901
RECYCLE  6   0.901   0.902
RECYCLE  7   0.902   0.902
RECYCLE  8   0.901   0.902
RECYCLE  9   0.901   0.902
Done2 with part2 (prediction)

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. rosettafoldna\_1.sh) for first step. For example:



```


#!/bin/bash
set -e
module load rosettafoldna
cd /data/$USER/rosettafoldna_test/
cp -r ${ROSETTAFOLDNA_TEST_DATA:-none}/* .
run_RF2NA_part1.sh test_o protein.fa R:RNA.fa


```

 Then create a batch input file (e.g. rosettafoldna\_2.sh) for the second step. For example:



```


#!/bin/bash
set -e
module load rosettafoldna
cd /data/$USER/rosettafoldna_test/
run_RF2NA_part2.sh test_o protein.fa R:RNA.fa


```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```

[user@biowulf]$ **sbatch --cpus-per-task=10 --mem=70g rosettafoldna\_1.sh**
1001
[user@biowulf]$ **sbatch --dependency=afterany:1001 --cpus-per-task=2 \
 --mem=10g --partition=gpu --gres=gpu:v100:1 rosettafoldna\_2.sh**
1002


```









