

document.querySelector('title').textContent = 'medaka on Biowulf';
medaka on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
 |


 medaka is a tool to create a consensus sequence from nanopore sequencing
data. This task is performed using neural networks applied from a pileup of
individual sequencing reads against a draft assembly. 



Documentation
* medaka on [GitHub](https://github.com/nanoporetech/medaka)
* medaka [Documentation](https://nanoporetech.github.io/medaka/)


Important Notes
* Module Name: medaka (see [the modules page](/apps/modules.html) for more information)
* medaka on biowulf can be run on CPU or GPU



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.


For this example we will use a P100 GPU.



```

[user@biowulf]$ **sinteractive --mem=32g -c12 --gres=lscratch:200,gpu:p100:1**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job
[user@cn3144]$ **cd /lscratch/$SLURM\_JOB\_ID**

```

As a first step, a draft assembly is created with pomoxis which includes a convenience
wrapper arount miniasm.



```

[user@cn3144]$ **module load pomoxis**
[user@cn3144]$ **cp -rL ${POMOXIS\_TEST\_DATA:-none} data**
[user@cn3144]$ **mini\_assemble -m 1 -i data/R9\_Ecoli\_K12\_MG1655\_lambda\_MinKNOW\_0.51.1.62.all.fasta \
 -o draft\_asm -p asm -t $SLURM\_CPUS\_PER\_TASK**
Copying FASTX input to workspace: data/R9_Ecoli_K12_MG1655_lambda_MinKNOW_0.51.1.62.all.fasta > draft_asm/asm.fa.gz
Skipped adapter trimming.
Skipped pre-assembly correction.
Overlapping reads...
[M::mm_idx_gen::26.495*2.10] collected minimizers
...etc...
Waiting for cleanup.
rm: cannot remove ‘shuffled*’: No such file or directory
rm: cannot remove ‘*paf*’: No such file or directory
Final assembly written to draft_asm/asm_final.fa. Have a nice day.

[user@cn3144]$ **ls -lh draft\_asm**
total 4.5M
-rw-r--r-- 1 user group 4.5M Sep 16 13:43 asm_final.fa
[user@cn3144]$ 

```

Polish the draft consensus with medaka. Note that a batch size (`-b`) of
150 works with our 16GB GPUs.


Note that there can be conflicts between pomoxis and medaka - they both include
their own mini\_align script - so pomoxis should be unloaded before running medaka.



```

[user@cn3144]$ **module load medaka/1.9.1**
[user@cn3144]$ **module unload pomoxis**
[user@cn3144]$ **medaka\_consensus -i data/R9\_Ecoli\_K12\_MG1655\_lambda\_MinKNOW\_0.51.1.62.all.fasta \
 -d draft\_asm/asm\_final.fa -o consensus -t $SLURM\_CPUS\_PER\_TASK -b 150**
This is medaka 1.5.0
Program    Version    Required   Pass
bcftools   1.14       1.11       True
bgzip      1.14       1.11       True
minimap2   2.17       2.11       True
samtools   1.14       1.11       True
tabix      1.14       1.11       True
Aligning basecalls to draft
...etc...
Running medaka stitch
[14:19:10 - DataIndex] Loaded sample-index from 1/1 (100.00%) of feature files.
[14:19:10 - Stitch] Processing utg000001l.
[14:19:17 - Stitch] Processing utg000002l.
[14:19:17 - Stitch] Processing utg000003l.
[14:19:17 - Stitch] Processing utg000004c.
Polished assembly written to consensus/consensus.fasta, have a nice day.


[user@cn3144]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. medaka.sh), which uses the input file 'medaka.in'. For example:



```

#!/bin/bash
module load medaka/1.9.1 || exit 1
medaka_consensus -i basecalls.fq -d draft_assembly.fa -o consensus -t $SLURM_CPUS_PER_TASK -b 150

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=12 --mem=24g medaka.sh
```







