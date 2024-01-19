

document.querySelector('title').textContent = 'pomoxis on Biowulf';
pomoxis on Biowulf


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



Pomoxis contains convenience wrappers around nanopore tools.



Documentation
* pomoxis on [GitHub](https://github.com/nanoporetech/pomoxis)
* pomoxis [Documentation](https://nanoporetech.github.io/pomoxis/index.html)


Important Notes
* Module Name: pomoxis (see [the modules page](/apps/modules.html) for more information)
* Some pomoxis tools are multithreaded. Please match the number of allocated CPUs and the number of threads
* Example files in `$POMOXIS_TEST_DATA`



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --gres=lscratch:50 --cpus-per-task=6 --mem=12g**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144]$ **cd /lscratch/$SLURM\_JOB\_ID**
[user@cn3144]$ **ml pomoxis**
[+] Loading minimap2, version 2.17...
[+] Loading miniasm 0.3.r179  ...
[+] Loading samtools 1.9  ...
[+] Loading racon 1.3.2  ...
[+] Loading seqkit  0.10.2
[+] Loading porechop  0.2.4
[+] Loading pomoxis  0.2.3

[user@cn3144]$ **cp -L $POMOXIS\_TEST\_DATA/\* .**
[user@cn3144]$ **ls -lh**
total 906M
-rw-r----- 1 user staff 4.5M Sep 12 09:12 NC_000913.3.fasta
-rw-r----- 1 user staff 901M Sep 12 09:12 R9_Ecoli_K12_MG1655_lambda_MinKNOW_0.51.1.62.all.fasta

```

The ecoli reads used in these examples were obtained from the 
[Loman Labs](http://lab.loman.net/2016/07/30/nanopore-r9-data-release/). The 
`mini_align` script uses minimap2 to align reads and runs some common post-alignment
processing steps like sorting and aligning the resulting bam file.



```


[user@cn3144]$ **mini\_align -r NC\_000913.3.fasta \
 -i R9\_Ecoli\_K12\_MG1655\_lambda\_MinKNOW\_0.51.1.62.all.fasta \
 -t $SLURM\_CPUS\_PER\_TASK -p ecoli**
Constructing minimap index.
[M::mm_idx_gen::0.178*1.01] collected minimizers
[M::mm_idx_gen::0.220*1.38] sorted minimizers
[M::main::0.274*1.31] loaded/built the index for 1 target sequence(s)
[M::mm_idx_stat] kmer size: 15; skip: 10; is_hpc: 0; #seq: 1
[M::mm_idx_stat::0.286*1.29] distinct minimizers: 838542 (98.18% are singletons); average occurrences: 1.034; average spacing: 5.352
[M::main] Version: 2.17-r941
[M::main] CMD: minimap2 -I 16G -x map-ont -d NC_000913.3.fasta.mmi NC_000913.3.fasta
[M::main] Real time: 0.292 sec; CPU: 0.376 sec; Peak RSS: 0.049 GB
[samfaipath] build FASTA index...
[M::main::0.070*1.02] loaded/built the index for 1 target sequence(s)
[M::mm_mapopt_update::0.086*1.01] mid_occ = 12
[M::mm_idx_stat] kmer size: 15; skip: 10; is_hpc: 0; #seq: 1
[M::mm_idx_stat::0.097*1.01] distinct minimizers: 838542 (98.18% are singletons); average occurrences: 1.034; average spacing: 5.352
[M::worker_pipeline::128.218*5.37] mapped 72602 sequences
[M::worker_pipeline::203.551*5.30] mapped 59810 sequences
[M::main] Version: 2.17-r941
[M::main] CMD: minimap2 -x map-ont -t 6 -a NC_000913.3.fasta.mmi R9_Ecoli_K12_MG1655_lambda_MinKNOW_0.51.1.62.all.fasta
[M::main] Real time: 203.564 sec; CPU: 1078.051 sec; Peak RSS: 3.832 GB
[bam_sort_core] merging from 0 files and 6 in-memory blocks...

```

De novo assembly with miniasm



```

[user@cn3144]$ **mini\_assemble -i R9\_Ecoli\_K12\_MG1655\_lambda\_MinKNOW\_0.51.1.62.all.fasta \
 -o ecoli\_assm -t $SLURM\_CPUS\_PER\_TASK -m 1 -c**
...much output...
[user@cn3144]$ **assess\_assembly -r NC\_000913.3.fasta \
 -i ecoli\_assm/reads\_final.fa -t $SLURM\_CPUS\_PER\_TASK**
...
  name     mean     q10      q50      q90
 err_ont  2.033%   1.783%   1.925%   2.332%
 err_bal  2.052%   1.797%   1.942%   2.357%
    iden  0.320%   0.257%   0.296%   0.344%
     del  0.821%   0.721%   0.788%   0.913%
     ins  0.913%   0.764%   0.858%   1.079%

#  Q Scores
  name     mean      q10      q50      q90
 err_ont  16.92    17.49    17.16    16.32
 err_bal  16.88    17.46    17.12    16.28
    iden  24.95    25.89    25.29    24.64
     del  20.85    21.42    21.03    20.40
     ins  20.40    21.17    20.67    19.67


[user@cn3144]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. pomoxis.sh), which uses the input file 'pomoxis.in'. For example:



```

#!/bin/bash
wd=$PWD
module load pomoxis/0.2.3 || exit 1
cd /lscratch/$SLURM_JOB_ID
cp -L $POMOXIS_TEST_DATA/R9_Ecoli_K12_MG1655_lambda_MinKNOW_0.51.1.62.all.fasta .
mini_assemble -i R9_Ecoli_K12_MG1655_lambda_MinKNOW_0.51.1.62.all.fasta \
    -o ecoli_assm -t $SLURM_CPUS_PER_TASK  -m 1 -c
mv ecoli_assm $wd

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=6 --gres=lscratch:50 --mem=10g pomoxis.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. pomoxis.swarm). For example:



```

mini_align -r NC_000913.3.fasta -i expt1.fastq -t $SLURM_CPUS_PER_TASK -p expt1
mini_align -r NC_000913.3.fasta -i expt2.fastq -t $SLURM_CPUS_PER_TASK -p expt2
mini_align -r NC_000913.3.fasta -i expt3.fastq -t $SLURM_CPUS_PER_TASK -p expt3

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f pomoxis.swarm -g 10 -t 6 --module pomoxis/0.2.3
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g #  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t #  Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module pomoxis  Loads the pomoxis module for each subjob in the swarm 
 | |
 | |
 | |








