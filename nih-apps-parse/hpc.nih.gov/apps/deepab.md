

document.querySelector('title').textContent = 'DeepAb on Biowulf';
DeepAb on Biowulf


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



DeepAb is antibody structure prediction using interpretable deep learning.



### References:


* JA Ruffolo, J Sulam, and JJ Gray 
 [**Antibody structure prediction using interpretable deep learning**](https://pubmed.ncbi.nlm.nih.gov/35199061/)
*Patterns (N Y). 2021 Dec 9;3(2):100406. doi: 10.1016/j.patter.2021.100406. eCollection 2022 Feb 11.*


Documentation
* [DeepAb Main Site](https://github.com/RosettaCommons/DeepAb)


Important Notes
* Module Name: deepab (see [the modules page](/apps/modules.html) for more information)
 * Multithreaded
 * GPU-dependent (only the predict.py script)
* Environment variables set 
	+ DEEPAB\_HOME* Example files in $DEEPAB\_HOME/data/sample\_files



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive --gres=gpu:p100:1 -c 8 --mem=20g**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load deepab**

[user@cn3144 ~]$ **predict.py --decoys 4 --num\_procs $SLURM\_CPUS\_ON\_NODE --use\_gpu --renumber \
--pred\_dir $(pwd)/preds --model\_dir $DEEPAB\_HOME/trained\_models/ensemble\_abresnet \
$DEEPAB\_HOME/data/sample\_files/4h0h.fasta**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. deepab.sh). For example:



```

#!/bin/bash
set -e
module load deepab
predict.py --decoys 4 --num_procs $SLURM_CPUS_ON_NODE --use_gpu --renumber --pred_dir $(pwd)/preds --model_dir $DEEPAB_HOME/trained_models/ensemble_abresnet $DEEPAB_HOME/data/sample_files/4h0h.fasta

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=8 --mem=20g --gres=gpu:p100:1 deepab.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. deepab.swarm). For example:



```

annotate_attention.py --renumber --cdr_loop h3 --model_file $DEEPAB_HOME/trained_models/ensemble_abresnet/rs0.pt --out_file output_1.pdb sample_1.pdb
annotate_attention.py --renumber --cdr_loop h3 --model_file $DEEPAB_HOME/trained_models/ensemble_abresnet/rs0.pt --out_file output_2.pdb sample_2.pdb
annotate_attention.py --renumber --cdr_loop h3 --model_file $DEEPAB_HOME/trained_models/ensemble_abresnet/rs0.pt --out_file output_3.pdb sample_3.pdb
annotate_attention.py --renumber --cdr_loop h3 --model_file $DEEPAB_HOME/trained_models/ensemble_abresnet/rs0.pt --out_file output_4.pdb sample_4.pdb

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f deepab.swarm -g 10g -t 2 --module deepab
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file).
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module deepab Loads the deepab module for each subjob in the swarm 
 | |
 | |
 | |








