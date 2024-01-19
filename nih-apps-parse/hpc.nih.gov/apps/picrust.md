

document.querySelector('title').textContent = 'PICRUSt on Biowulf';
PICRUSt on Biowulf


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



PICRUSt (pronounced “pie crust”) is a bioinformatics software package designed to predict metagenome functional content from marker gene (e.g., 16S rRNA) surveys and full genomes.



### References:


* [Douglas, Gavin M., et al. "PICRUSt2 for prediction of metagenome functions." *Nature Biotechnology* 38.6 (2020): 685-688.](https://www.nature.com/articles/s41587-020-0548-6)


Documentation
* [PICRUSt Main Site](https://huttenhower.sph.harvard.edu/picrust)
* [PICRUSt on GitHub](https://github.com/picrust/picrust2)


Important Notes
* Module Name: picrust (see [the modules page](/apps/modules.html) for more information)
 * Remember to set the processes option/argument pair like so: -p $SLURM\_CPUS\_PER\_TASK* Environment variables set 
	+ PICRUST\_HOME* Example files in $PICRUST\_HOME/tutorials* Recent versions of PICRUSt on Biowulf are technically PICRUSt2



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf ~]$ **sinteractive -c8 --mem=16g --gres=lscratch:10**
salloc: Pending job allocation 30152043
salloc: job 30152043 queued and waiting for resources
salloc: job 30152043 has been allocated resources
salloc: Granted job allocation 30152043
salloc: Waiting for resource configuration
salloc: Nodes cn0890 are ready for job
srun: error: x11: no local DISPLAY defined, skipping
error: unable to open file /tmp/slurm-spank-x11.30152043.0
slurmstepd: error: x11: unable to read DISPLAY value

[user@cn0890 ~]$ **cd /lscratch/$SLURM\_JOB\_ID**

[user@cn0890 30152043]$ **module load picrust/2.4.2**
[+] Loading picrust  2.4.2  on cn0890
[+] Loading singularity  3.8.5-1  on cn0890

[user@cn0890 30152043]$ **cp $PICRUST\_HOME/tutorials/chemerin\_16S.zip .**

[user@cn0890 30152043]$ **unzip chemerin\_16S.zip**
Archive:  chemerin_16S.zip
  inflating: chemerin_16S/metadata.tsv
  inflating: chemerin_16S/seqs.fna
  inflating: chemerin_16S/table.biom

[user@cn0890 30152043]$ **picrust2\_pipeline.py -s chemerin\_16S/seqs.fna -i chemerin\_16S/table.biom -o \
 picrust2\_out\_pipeline -p $SLURM\_CPUS\_PER\_TASK**




All ASVs were below the max NSTI cut-off of 2.0 and so all were retained for downstream analyses.

All ASVs were below the max NSTI cut-off of 2.0 and so all were retained for downstream analyses.



[user@cn0890 30152043]$ **ls picrust2\_out\_pipeline/**
EC_metagenome_out    intermediate       KO_predicted.tsv.gz               out.tre
EC_predicted.tsv.gz  KO_metagenome_out  marker_predicted_and_nsti.tsv.gz  pathways_out

[user@cn0890 30152043]$ **exit**
exit
salloc: Relinquishing job allocation 30152043

[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. picrust.sh). For example:



```

#!/bin/bash
set -e
module load picrust
picrust2_pipeline.py -s seqs.fna -i table.biom -o output -p $SLURM_CPUS_PER_TASK

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] picrust.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. picrust.swarm). For example:



```

picrust2_pipeline.py -s seqs1.fna -i table1.biom -o output -p $SLURM_CPUS_PER_TASK
picrust2_pipeline.py -s seqs2.fna -i table2.biom -o output -p $SLURM_CPUS_PER_TASK
picrust2_pipeline.py -s seqs3.fna -i table3.biom -o output -p $SLURM_CPUS_PER_TASK
picrust2_pipeline.py -s seqs4.fna -i table4.biom -o output -p $SLURM_CPUS_PER_TASK

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f picrust.swarm [-g #] [-t #] --module picrust
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module picrust Loads the picrust module for each subjob in the swarm 
 | |
 | |
 | |








