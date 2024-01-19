

document.querySelector('title').textContent = 'colabfold on Biowulf';

 .hl { background-color: #ffff99; }
 code {background-color: #eeeeee; border: 1px solid #bbbbbb; padding: 0;}
 dt {font-weight: bold; margin-top: 5px;}
 dd {padding-left: 2px; border-left: 1px solid #bbbbbb;}
 .btt {border: 1px solid silver;
 background-color: white;
 padding: 5px;
 position: relative;
 margin: 5px 0px 10px 10px;
 float: right;
 top: -25px;
 left: 10px;
 }

colabfold on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
 |



This module provides the batch scripts of the ColabFold implementation of alphafold.



### References:


* Mirdita M, Schütze K, Moriwaki Y, Heo L, Ovchinnikov S and Steinegger M. 
 *ColabFold: Making protein folding accessible to all.*
 Nature Methods (2022) doi: 10.1038/s41592-022-01488-1. 
 [PubMed](https://pubmed.ncbi.nlm.nih.gov/35637307/) | 
 [PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9184281/) | 
 [Journal](https://www.nature.com/articles/s41592-022-01488-1)


Documentation
* colabfold on [GitHub](https://github.com/sokrypton/ColabFold)


Important Notes
* Module Name: colabfold (see [the modules page](/apps/modules.html) for more information)
* The MSA generation is a multithreaded process. Model generation uses GPU
* Example files in `$COLABFOLD_TEST_DATA`
* Reference data in `/fdb/colabfold/`
* Currently including templates is not yet supported



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) for
the generation of the multiple sequence applications (MSAs). This session
requires at least 128GB of memory. Note that
*colabfold\_search is optimized for running many query sequences* in a
single job. Never run single sequences in the search part of
the analysis. It is inefficient and if you run many will strain the file
system. For example: Below we will create alignments for 17 polymerase
related proteins of the mimivirus genome which takes about 100 minutes.
Creating MSAs for all 979 proteins of the mimivirus genome take takes about 285
minutes or only 2.8x longer for 50x more proteins.

Note that colabfold\_search
treats each protein in a fasta file as a separate monomer.



```

[user@biowulf]$ **sinteractive --mem=128G --cpus-per-task=16 --gres=lscratch:100**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144]$ **module load colabfold**
[user@cn3144]$ **cp $COLABFOLD\_TEST\_DATA/mimi\_poly.fa .**
[user@cn3144]$ **colabfold\_search --threads $SLURM\_CPUS\_PER\_TASK \
 mimi\_poly.fa $COLABFOLD\_DB mimi\_poly**
[...much output...]
[user@cn3144]$ **ls -lh mimi\_poly**
total 66M
-rw-r--r-- 1 user group  17M Sep 22 16:14 0.a3m
-rw-r--r-- 1 user group 2.7M Sep 22 16:14 10.a3m
-rw-r--r-- 1 user group  92K Sep 22 16:14 11.a3m
-rw-r--r-- 1 user group 296K Sep 22 16:14 12.a3m
-rw-r--r-- 1 user group  33K Sep 22 16:14 13.a3m
-rw-r--r-- 1 user group 6.4M Sep 22 16:14 14.a3m
-rw-r--r-- 1 user group  42K Sep 22 16:14 15.a3m
-rw-r--r-- 1 user group 7.1K Sep 22 16:14 16.a3m
-rw-r--r-- 1 user group  23M Sep 22 16:14 1.a3m
-rw-r--r-- 1 user group  12M Sep 22 16:14 2.a3m
-rw-r--r-- 1 user group 797K Sep 22 16:14 3.a3m
-rw-r--r-- 1 user group 200K Sep 22 16:14 4.a3m
-rw-r--r-- 1 user group 580K Sep 22 16:14 5.a3m
-rw-r--r-- 1 user group 273K Sep 22 16:14 6.a3m
-rw-r--r-- 1 user group 160K Sep 22 16:14 7.a3m
-rw-r--r-- 1 user group 464K Sep 22 16:14 8.a3m
-rw-r--r-- 1 user group 2.6M Sep 22 16:14 9.a3m
[user@cn3144]$ **exit**

```

Now we can build the structure predictions for the 17 proteins above. This step requires
GPU. Each of the 17 proteins in this example is modeled as a monomer. See below for
multimer predictions. With the default settings on a A100 model generation for the 17
proteins takes approximately 4h



```

[user@biowulf]$ **sinteractive --mem=48G --cpus-per-task=8 --gres=lscratch:100,gpu:a100:1**
salloc.exe: Pending job allocation 46116227
salloc.exe: job 46116227 queued and waiting for resources
salloc.exe: job 46116227 has been allocated resources
salloc.exe: Granted job allocation 46116227
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job
[user@cn2144]$ **module load colabfold**
[user@cn2144]$ **colabfold\_batch --amber --use-gpu-relax \
 mimi\_poly mimi\_poly\_models**
2022-09-26 09:54:12,026 Running colabfold 1.3.0 (22671664ac2c9dcb30086c3e654414d950ccb297)
2022-09-26 09:54:19,458 Found 6 citations for tools or databases
2022-09-26 09:54:24,671 Query 1/17: 12 (length 73)
2022-09-26 09:54:24,720 Running model_3
2022-09-26 09:54:50,065 model_3 took 23.3s (3 recycles) with pLDDT 86.1 and ptmscore 0.622
2022-09-26 09:55:15,397 Relaxation took 15.9s
2022-09-26 09:55:15,400 Running model_4
2022-09-26 09:55:24,884 model_4 took 8.6s (3 recycles) with pLDDT 89 and ptmscore 0.668
2022-09-26 09:55:31,969 Relaxation took 5.7s
[...snip...]
[user@cn2144]$ **ls -lh mimi\_poly\_models | head**
total 570M
-rw-r--r-- 1 user group  17M Sep 26 15:36 0.a3m
-rw-r--r-- 1 user group 142K Sep 26 16:19 0_coverage.png
-rw-r--r-- 1 user group    0 Sep 26 16:19 0.done.txt
-rw-r--r-- 1 user group 743K Sep 26 16:19 0_PAE.png
-rw-r--r-- 1 user group 196K Sep 26 16:19 0_plddt.png
-rw-r--r-- 1 user group  18M Sep 26 16:19 0_predicted_aligned_error_v1.json
-rw-r--r-- 1 user group 1.5M Sep 26 16:19 0_relaxed_rank_1_model_3.pdb
-rw-r--r-- 1 user group 1.5M Sep 26 16:19 0_relaxed_rank_2_model_4.pdb
-rw-r--r-- 1 user group 1.5M Sep 26 16:19 0_relaxed_rank_3_model_2.pdb
[user@cn2144]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf]$

```

For each protein the output includes some diagnostic plots, the predicted alignment error
of the ranked models in json, as well as the ranked relaxed and unrelaxed models in PDB format.
As an example below see the predicted best ranked model (rank 1) for YP\_003986740 (DNA dependent
RNA polymerase subunit B) along with some of the diagnostic plots.




![colabfold predicted structure for YP_003986740](/images/colabfold_fig1.png)
**Figure 1**. **(A)** The highest confidence
 structure prediction (`0_relaxed_rank_1_model_3.pdb`) of protein
 YP\_003986740 colored by pLDDT. **(B)** pLDDT for all 5 models.
 **(C)** Coverage of YP\_003986740 in the MSA generated by colabfold\_search
 with mmseqs2. **(D)** Predicted aligned error (PAE) for all 5 models.


Commandline for generating two predictions for each of the 5 models per protein but only
relax the top 5 results



```

[user@cn2144]$ **colabfold\_batch --amber --use-gpu-relax --num-relax 5 --num-seeds 2 ...**

```

#### Predicting multimers with colabfold batch tools


Multimers are predicted by concatenating all sequences of a multimer separated by ':' into a
single sequence



Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. colabfold-search.sh) to create the MSAs:



```

#!/bin/bash
module load colabfold
cp $COLABFOLD_TEST_DATA/mimi_poly.fa .
colabfold_search --threads $SLURM_CPUS_PER_TASK \
    mimi_poly.fa $COLABFOLD_DB mimi_poly

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=16 --mem=128g --gres=lscratch:100 colabfold-search.sh
```

Then build the corresponding models



```

#!/bin/bash
module load colabfold
colabfold_batch --amber --use-gpu-relax mimi_poly mimi_poly_models

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=8 --mem=48g --gres=lscratch:100,gpu:a100:1 --partition=gpu colabfold-batch.sh
```







