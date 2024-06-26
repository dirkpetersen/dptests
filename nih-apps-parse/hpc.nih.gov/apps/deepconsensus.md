

document.querySelector('title').textContent = 'DeepConsensus on Biowulf';
DeepConsensus on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Benchmarks](#metrics)
[Interactive job](#int) 
[Batch job](#sbatch) 
[Swarm of jobs](#swarm) 
 |



DeepConsensus is used to correct PacBio sequence data. It helps convert circular consensus sequence (CCS) data into corrected FASTQ files for downstream analysis. According to the authors:




> 
>  DeepConsensus ... uses a unique alignment-based loss to train a gap-aware transformer-encoder (GATE) for sequence correction. Compared to pbccs, DeepConsensus reduces read errors in the same dataset by 42%. 
> 


### References:


* Gunjan Baid, et al. 
 [**DeepConsensus: Gap-Aware Sequence Transformers for Sequence Correction**](https://doi.org/10.1101/2021.08.31.458403)
*bioRxiv 2021.08.31.458403.*


Documentation
* [deepconsensus GitHub](https://github.com/google/deepconsensus)
* [Quick Start v0.3](https://github.com/google/deepconsensus/blob/r0.3/docs/quick_start.md)


Important Notes
* Module Name: deepconsensus (see [the modules page](/apps/modules.html) for more information)
* Multithreaded, GPU capable
* Environment variables set:
+ DEEPCONSENSUS\_EXAMPLE

* Example files in /usr/local/apps/deepconsensus/examples


GPU Performance Benchmarks
The following benchmarks were generated by running the quick start v0.2 example ([as posted on GitHub](https://github.com/google/deepconsensus/blob/r0.2/docs/quick_start.md)) across our GPUs using the deepconsensus/0.2.0-gpu module. The jobs allocated 1 GPU, 16 CPUs, and 16G memory. The --batch\_zmws option was set to 100.


Listed below are the total times in seconds to process 1000 ZMWs as reported by deepconsensus. Compare these values with those posted by the developers [here](https://github.com/google/deepconsensus/blob/main/docs/runtime_metrics.md). Note that our installations of deepconsensus are container-based.


* k20x: 790s
* k80: 696s
* p100: 405s
* v100: 351s
* v100x: 343s
* a100: 164s


Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive --gres=lscratch:10,gpu:p100:1 --mem=16G --cpus-per-task=16**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load deepconsensus/0.3.1-gpu**

[user@cn3144 ~]$ **cd /lscratch/$SLURM\_JOB\_ID**

[user@cn3144 ~]$ **cp -r $DEEPCONSENSUS\_EXAMPLE .**

[user@cn3144 ~]$ **cd examples**

[user@cn3144 ~]$ **shard\_id=00001-of-00002**

[user@cn3144 ~]$  **pbindex n1000.subreads.bam**

[user@cn3144 ~]$ **ccs --min-rq=0.88 \
 -j 16 \
 --chunk=1/2 \
 n1000.subreads.bam \
 "${shard\_id}.ccs.bam"**

[user@cn3144 ~]$ **actc -j "$(nproc)" \
 n1000.subreads.bam \
 "${shard\_id}.ccs.bam" \
 "${shard\_id}.subreads\_to\_ccs.bam"**

[user@cn3144 ~]$ **deepconsensus run \
 --cpus 16 \
 --subreads\_to\_ccs=${shard\_id}.subreads\_to\_ccs.bam \
 --ccs\_bam=${shard\_id}.ccs.bam \
 --checkpoint=model/checkpoint \
 --output=${shard\_id}.output.fastq**
I0908 23:02:15.754746 46912500000576 quick_inference.py:727] Using multiprocessing: cpus is 16.
I0908 23:02:15.755900 46912500000576 quick_inference.py:459] Loading model/checkpoint
2022-09-08 23:02:15.776381: I tensorflow/core/platform/cpu_feature_guard.cc:193] This TensorFlow binary is optimized with oneAPI Deep Neural Network Library (oneDNN) to use the following CPU instructions in performance-critical operations:  AVX2 FMA
To enable them in other operations, rebuild TensorFlow with the appropriate compiler flags.
2022-09-08 23:02:18.999404: I tensorflow/core/common_runtime/gpu/gpu_device.cc:1532] Created device /job:localhost/replica:0/task:0/device:GPU:0 with 15401 MB memory:  -> device: 0, name: Tesla P100-PCIE-16GB, pci bus id: 0000:0d:00.0, compute capability: 6.0
I0908 23:02:19.214797 46912500000576 networks.py:358] Condensing input.
Model: "encoder_only_learned_values_transformer"
...
I0908 23:02:22.412791 46912500000576 model_utils.py:178] 1 GPUs being used.
I0908 23:02:22.413068 46912500000576 model_utils.py:179] Per-replica batch-size is 32768.
I0908 23:02:22.413206 46912500000576 model_utils.py:181] Global batch-size is 32768.
I0908 23:02:22.413361 46912500000576 model_utils.py:231] Setting hidden size to transformer_input_size.
I0908 23:02:22.413509 46912500000576 quick_inference.py:484] Finished initialize_model.
I0908 23:02:22.416144 46912500000576 quick_inference.py:738] Model setup took 6.661123037338257 seconds.
I0908 23:03:43.689472 46912500000576 quick_inference.py:623] Example summary: ran model=21414 (66.91%; 45.476s) skip=10591 (33.09%; 2.068s) total=32005.
I0908 23:03:44.986394 46912500000576 quick_inference.py:662] Processed a batch of 100 ZMWs in 66.135 seconds
I0908 23:03:45.069340 46912500000576 quick_inference.py:779] Processed 100 ZMWs in 82.653 seconds
I0908 23:04:50.641679 46912500000576 quick_inference.py:623] Example summary: ran model=18669 (71.82%; 38.282s) skip=7326 (28.18%; 1.532s) total=25995.
I0908 23:04:51.698215 46912500000576 quick_inference.py:662] Processed a batch of 78 ZMWs in 53.419 seconds
I0908 23:04:51.901666 46912500000576 quick_inference.py:798] Processed 178 ZMWs in 149.485 seconds
I0908 23:04:51.901839 46912500000576 quick_inference.py:800] Outcome counts: OutcomeCounter(empty_sequence=0, only_gaps_and_padding=0, failed_quality_filter=3, failed_length_filter=0, success=175)

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. deepconsensus.sh). For example:



```

#!/bin/bash
set -e
module load deepconsensus/0.3.1-gpu
cd /lscratch/$SLURM_JOB_ID
ccs --all -j $SLURM_CPUS_PER_TASK /data/$USER/subreads.bam ccs.bam
actc -j $SLURM_CPUS_PER_TASK /data/$USER/subreads.bam ccs.bam subreads_to_ccs.bam
deepconsensus run \
              --subreads_to_ccs=subreads_to_ccs.bam \
              --ccs_fasta=ccs.fasta \
              --output=output.fastq \
              --batch_zmws=100 \
              --cpus $SLURM_CPUS_PER_TASK
mv output.fastq /data/$USER/

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#G] [--gres=lscratch:#,gpu:k80:1] deepconsensus.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. deepconsensus.swarm). For example:



```

deepconsensus run --subreads_to_ccs= ... --output=output_1.fastq --batch_zmws=100 --cpus $SLURM_CPUS_PER_TASK
deepconsensus run --subreads_to_ccs= ... --output=output_2.fastq --batch_zmws=100 --cpus $SLURM_CPUS_PER_TASK
deepconsensus run --subreads_to_ccs= ... --output=output_3.fastq --batch_zmws=100 --cpus $SLURM_CPUS_PER_TASK
deepconsensus run --subreads_to_ccs= ... --output=output_4.fastq --batch_zmws=100 --cpus $SLURM_CPUS_PER_TASK

```

Make sure to replace ... with other required deepconsensus options and flags. Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f deepconsensus.swarm [-g #] [-t #] --gres=gpu:k20x:1 --module deepconsensus/0.2.0-gpu
```

where


|  |  |
| --- | --- |
| -g *#*  | Number of Gigabytes of memory required for each process (1 line in the swarm command file) |
| -t *#* | Number of threads/CPUs required for each process (1 line in the swarm command file). |
| --gre=gpu:k20x:1 | Allocates 1 k20x GPU for each process (1 line in the swarm command file). |
| --module deepconsensus/0.2.0-gpu | Loads the deepconsensus/0.2.0-gpu module for each subjob in the swarm |






