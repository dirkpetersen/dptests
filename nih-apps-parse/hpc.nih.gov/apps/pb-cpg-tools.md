

document.querySelector('title').textContent = 'pb-CpG-tools on Biowulf';
pb-CpG-tools on Biowulf


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



A collection of tools for the analysis of CpG/5mC data from PacBio HiFi reads aligned to a reference genome (e.g., an aligned BAM). To use these tools, the HiFi reads should already contain 5mC base modification tags, generated on-instrument or by using primrose. The aligned BAM should also be sorted and indexed.



Documentation
* [pb-CpG-tools Github](https://github.com/PacificBiosciences/pb-CpG-tools)


Important Notes
* Module Name: pb-cpg-tools (see [the modules page](/apps/modules.html) for more information)
* Multithreaded
* Environment variables set 
	+ CPG\_TEST\_DATA
	+ CPG\_PILEUP\_MODEL* Test data files in $CPG\_TEST\_DATA
* The models distributed with the git repository can be copied from $CPG\_PILEUP\_MODEL directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive --mem=32G -c8 --gres=lscratch:20**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load pb-cpg-tools**

[user@cn3144 ~]$ **cd /lscratch/$SLURM\_JOB\_ID**

[user@cn3144 ~]$ **cp $CPG\_TEST\_DATA/\*bam\* .**

[user@cn3144 ~]$ **aligned\_bam\_to\_cpg\_scores \
 --bam HG002.GRCh38.haplotagged.truncated.bam \
 --output-prefix test \
 --model $CPG\_PILEUP\_MODEL/pileup\_calling\_model.v1.tflite \
 --threads $SLURM\_CPUS\_PER\_TASK**

[2023-07-24][14:31:38][aligned_bam_to_cpg_scores][INFO] Starting aligned_bam_to_cpg_scores
[2023-07-24][14:31:38][aligned_bam_to_cpg_scores][INFO] cmdline: aligned_bam_to_cpg_scores --bam HG002.GRCh38.haplotagged.truncated.bam --output-prefix test --model /usr/local/apps/pb-cpg-tools/2.3.1/models/pileup_calling_model.v1.tflite --threads 8
[2023-07-24][14:31:38][aligned_bam_to_cpg_scores][INFO] Running on 8 threads
[2023-07-24][14:31:38][aligned_bam_to_cpg_scores][INFO] Processing alignment file 'HG002.GRCh38.haplotagged.truncated.bam'
[2023-07-24][14:32:16][aligned_bam_to_cpg_scores][INFO] Finished processing alignment files.
[2023-07-24][14:32:16][aligned_bam_to_cpg_scores][INFO] Writing hap2 site methylation to bed file: 'test.hap2.bed'
[2023-07-24][14:32:16][aligned_bam_to_cpg_scores][INFO] Writing hap1 site methylation to bed file: 'test.hap1.bed'
[2023-07-24][14:32:16][aligned_bam_to_cpg_scores][INFO] Writing combined site methylation to bed file: 'test.combined.bed'
[2023-07-24][14:32:16][aligned_bam_to_cpg_scores][INFO] Writing combined site methylation to bigwig file: 'test.combined.bw'
[2023-07-24][14:32:16][aligned_bam_to_cpg_scores][INFO] Writing hap2 site methylation to bigwig file: 'test.hap2.bw'
[2023-07-24][14:32:16][aligned_bam_to_cpg_scores][INFO] Writing hap1 site methylation to bigwig file: 'test.hap1.bw'
[2023-07-24][14:32:17][aligned_bam_to_cpg_scores][INFO] aligned_bam_to_cpg_scores completed. Total Runtime: 00:00:39.426

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. pb-cpg-tools.sh). For example:



```

#!/bin/bash
set -e
module load pb-cpg-tools
cd /lscratch/$SLURM_JOB_ID
cp $CPG_TEST_DATA/*bam* .
aligned_bam_to_cpg_scores \
            --bam HG002.GRCh38.haplotagged.truncated.bam \
            --output-prefix test \
            --model $CPG_PILEUP_MODEL/pileup_calling_model.v1.tflite \
            --threads $SLURM_CPUS_PER_TASK

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --mem=32g --cpus-per-task=8 --gres=lscratch:20  pb-cpg-tools.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. pb-cpg-tools.swarm). For example:



```

aligned_bam_to_cpg_scores --bam input1.bam --output-prefix out1 --model v1.tflite --threads $SLURM_CPUS_PER_TASK
aligned_bam_to_cpg_scores --bam input2.bam --output-prefix out2 --model v1.tflite --threads $SLURM_CPUS_PER_TASK
aligned_bam_to_cpg_scores --bam input3.bam --output-prefix out3 --model v1.tflite --threads $SLURM_CPUS_PER_TASK
aligned_bam_to_cpg_scores --bam input4.bam --output-prefix out4 --model v1.tflite --threads $SLURM_CPUS_PER_TASK

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f pb-cpg-tools.swarm [-g #] [-t #] [--gres=lscratch:#] --module pb-cpg-tools
```

where


|  |  |
| --- | --- |
| -g *#* | Number of Gigabytes of memory required for each process (1 line in the swarm command file) |
| -t *#* | Number of threads/CPUs required for each process (1 line in the swarm command file). |
| --gres=lscratch:# | lscratch amount in GB allocated for each process (1 line in the swarm command file). |
| --module pb-cpg-tools | Loads the pb-cpg-tools module for each subjob in the swarm |








