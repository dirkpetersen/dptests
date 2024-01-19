

document.querySelector('title').textContent = 'megalodon on Biowulf';
megalodon on Biowulf


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



Megalodon is a research command line tool to extract high accuracy modified base and sequence variant calls from raw nanopore reads by anchoring the information rich basecalling neural network output to a reference genome/transcriptome.



### References:


* [Gouil, Quentin, and Andrew Keniry. "Latest techniques to study DNA methylation." *Essays in Biochemistry* 63.6 (2019): 639-648.](https://portlandpress.com/essaysbiochem/article/63/6/639/221208/Latest-techniques-to-study-DNA-methylation)


Documentation
* [megalodon on GitHub](https://github.com/nanoporetech/megalodon)
* [megalodon Documentation](https://nanoporetech.github.io/megalodon/)


Important Notes
* Module Name: megalodon (see [the modules page](/apps/modules.html) for more information)
 * This is a GPU-only application and it will only run on Pascal or higher architecture. It will NOT run on Kepler (k20 or k80) GPUs. 
 * Example files in /usr/local/apps/megalodon/<ver>/example* Use --guppy-server-path ${MEGALODON\_GUPPY\_PATH} in your megalodon call. It points to the correct executable *within the container*.



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. The data in this example do not produce useful results but are for illustrative purposes only.   
Sample session (user input in **bold**):



```

[user@biowulf ~]$ **sinteractive -c12 --mem=32g --gres=lscratch:10,gpu:p100:1**
salloc.exe: Pending job allocation 5091510
salloc.exe: job 5091510 queued and waiting for resources
salloc.exe: job 5091510 has been allocated resources
salloc.exe: Granted job allocation 5091510
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn2369 are ready for job
srun: error: x11: no local DISPLAY defined, skipping

[user@cn2369 ~]$ **cd /lscratch/$SLURM\_JOB\_ID**

[user@cn2369 5091510]$ **module load megalodon**
[+] Loading megalodon  2.5.0  on cn2369
[+] Loading singularity  3.8.1-1  on cn2369

[user@cn2369 5091510]$ **cp -r /usr/local/apps/megalodon/2.5.0/TEST\_DATA/ .**

[user@cn2369 5091510]$ **cd TEST\_DATA/**

[user@cn2369 example]$ **megalodon fast5/ \
 --guppy-server-path ${MEGALODON\_GUPPY\_PATH} \
 --outputs basecalls mappings mod\_mappings mods \
 --output-directory mega\_out \
 --reference /fdb/igenomes/Homo\_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa \
 --devices ${CUDA\_VISIBLE\_DEVICES} \
 --processes ${SLURM\_CPUS\_PER\_TASK} \
 --overwrite \
 --guppy-config res\_dna\_r941\_min\_modbases\_5mC\_CpG\_v001.cfg**
******************** WARNING: "mods" output requested, so "per_read_mods" will be added to outputs. ********************
[12:27:45] Loading guppy basecalling backend
[2022-07-12 12:27:49.969485] [0x00002aaaf57a0700] [info]    Connecting to server as ''
[2022-07-12 12:27:49.972970] [0x00002aaaf57a0700] [info]    Connected to server as ''. Connection id: e96656bc-1a2a-4e39-8ea3-0dccc1f0757a
[12:27:50] Loading reference
[12:29:32] Loaded model calls canonical alphabet ACGT and modified bases m=5mC (alt to C)
[12:29:33] Preparing workers to process reads
[12:29:34] Processing reads
Full output or empty input queues indicate I/O bottleneck
3 most common unsuccessful processing stages:
    -----
[2022-07-12 12:29:34.642849] [0x00002aaaed59f700] [info]    Connecting to server as ''
[2022-07-12 12:29:34.646135] [0x00002aaaed59f700] [info]    Connected to server as ''. Connection id: 4b6bae96-7baa-44a2-84fa-a8ab743d4efa, ?reads/s]
...
[2022-07-12 12:29:37.269726] [0x00002aaaed59f700] [info]    Connecting to server as ''
[2022-07-12 12:29:37.286846] [0x00002aaaed59f700] [info]    Connected to server as ''. Connection id: 03633509-8c2e-45a1-b4a9-feda7170a5ca
    99.6% (   7871 reads) : No alignment
    -----
    -----
Read Processing: 100%|██████████████████████████████| 8000/8000 [02:03<00:00, 64.97reads/s, samples/s=2.4e+6]
 input queue capacity extract_signal      :   0%|                                                   | 0/10000
output queue capacity basecalls           :   0%|                                                   | 0/10000
output queue capacity mappings            :   0%|                                                   | 0/10000
output queue capacity per_read_mods       :   0%|                                                   | 0/10000
[12:31:37] Unsuccessful processing types:
    99.7% (   7972 reads) : No alignment
[12:31:38] Spawning modified base aggregation processes
[12:31:40] Aggregating 5946 per-read modified base statistics
[12:31:40] NOTE: If this step is very slow, ensure the output directory is located on a fast read disk (e.g. local SSD). Aggregation can be restarted using the `megalodon_extras aggregate run` command
Mods: 100%|███████████████████████████████████████████████| 5946/5946 [00:00<00:00, 27328.23 per-read calls/s][12:31:41] Mega Done

[user@cn2369 example]$ **exit**
exit
salloc.exe: Relinquishing job allocation 5091510

[user@biowulf ~]$ 

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. megalodon.sh). For example:



```

#!/bin/bash
set -e
module load megalodon
megalodon fast5/ \
     --guppy-server-path ${MEGALODON_GUPPY_PATH} \
     --outputs basecalls mappings mod_mappings mods \
     --output-directory mega_out \
     --reference /fdb/igenomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa \
     --devices ${CUDA_VISIBLE_DEVICES} \
     --processes ${SLURM_CPUS_PER_TASK} \
     --overwrite

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] [--gres=cpu:#] megalodon.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. megalodon.swarm). For example:



```

megalodon dir1/ \
     --guppy-server-path ${MEGALODON_GUPPY_PATH} \
     --outputs basecalls mappings mod_mappings mods \
     --output-directory out1 \
     --reference /fdb/igenomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa \
     --devices ${CUDA_VISIBLE_DEVICES} \
     --processes ${SLURM_CPUS_PER_TASK}
megalodon dir2/ \
     --guppy-server-path ${MEGALODON_GUPPY_PATH} \
     --outputs basecalls mappings mod_mappings mods \
     --output-directory out2 \
     --reference /fdb/igenomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa \
     --devices ${CUDA_VISIBLE_DEVICES} \
     --processes ${SLURM_CPUS_PER_TASK}
megalodon dir3/ \
     --guppy-server-path ${MEGALODON_GUPPY_PATH} \
     --outputs basecalls mappings mod_mappings mods \
     --output-directory out3 \
     --reference /fdb/igenomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa \
     --devices ${CUDA_VISIBLE_DEVICES} \
     --processes ${SLURM_CPUS_PER_TASK}

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f megalodon.swarm [-g #] [-t #] [--gres <gpu>:#] --module megalodon
```

where


|  |  |  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --gres *#* GPU type and number to use for each subjob.
 | --module megalodon Loads the megalodon module for each subjob in the swarm 
 | |
 | |
 | |
 | |








