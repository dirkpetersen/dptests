

document.querySelector('title').textContent = "alphapulldown";

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

alphapulldown on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
 |


 Streamlines protein-protein interaction screens and high-throughput
modelling of higher-order oligomers using AlphaFold-Multimer


### References:


* D. Yu, G. Chojnowski, M. Rosenthal, and J. Josinski. 
 *AlphaPulldown—a python package for protein–protein interaction screens using AlphaFold-Multimer*
 Bioinformatics 2023.
 [PubMed](https://pubmed.ncbi.nlm.nih.gov/36413069/) | 
 [PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/pmid/36413069/) | 
 [Journal](https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/btac749)


Documentation
* AlphaPulldown on [GitHub](https://github.com/KosinskiLab/AlphaPulldown)


Important Notes
* Module Name: `alphapulldown` (see [the modules page](/apps/modules.html) 
 for more information)
* Some steps in an analysis require GPU
* Parallelizing the default alignment generation method adopted from alphafold has to be done with
 care as it can easily overload file systems and generate more alignments than can be analyzed with
 available GPUs in a reasonable amount of time. Limit parallelization to 10 concurrent jobs for this step or
 use the mmseqs based `colabfold_search` from the
 [colabfold](https://hpc.nih.gov/apps/colabfold.html) module.
* Example files in `$ALPHAPULLDOWN_TEST_DATA`



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
In the following example we will screen a bait against a set of candidate
interaction partners (i.e. *pulldown* mode).


### Generating multiple sequence alignments (MSAs)


Alignment generation does not require GPUs and can be done with the either
the alignment pipeline from alphafold (included in this module) or colabfold
(separate module). The mmseqs based colabfold method is considerably faster.
Allocate an [interactive session](/docs/userguide.html#int) without
GPUs for the alignment step.


#### Alphafold based alignment pipeline



```

[user@biowulf]$ **sinteractive --mem=100g --cpus-per-task=8 --gres=lscratch:100 --time=1-12**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144]$ **module load alphapulldown/0.30.7**
[user@cn3144]$ **cp ${ALPHAPULLDOWN\_TEST\_DATA:-none}/\* .**
[user@cn3144]$ **ls -lh**
-rw-r--r-- 1 user group 1.1K Sep 12 10:50 bait.fasta
-rw-r--r-- 1 user group    7 Sep 12 20:25 bait.txt
-rw-r--r-- 1 user group  145 Sep 12 20:25 candidates.txt
-rw-r--r-- 1 user group 7.4K Sep 12 10:49 candidates.fasta
[user@cn3144]$ ### alphafold based alignment pipeline - not parallelized
[user@cn3144]$ **create\_individual\_features.py \
 --fasta\_paths=bait.fasta,candidates.fasta \
 --save\_msa\_files=True \
 --output\_dir=pulldown \
 --use\_precomputed\_msas=False \
 --max\_template\_date=2023-01-01 \
 --skip\_existing=True**
[...much output...]
[user@cn3144]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf]$

```

MSA generation takes 35h in the example above. Note that in an interactive
job it cannot be parallelized. To parallelize across subunits please see the
section on running a [Batch job](#sbatch) . Note that not saving msa
files saves space (12.7GB/15GB in this example) but precludes their re-use in
the future.


#### Colabfold based alignment pipeline



```

[user@biowulf]$ **sinteractive --mem=128g --cpus-per-task=16 --gres=lscratch:100**
[user@cn3144]$ **cat bait.fasta candidates.fasta > for\_colabfold.fa**
[user@cn3144]$ **module load colabfold alphapulldown**
[user@cn3144]$ **colabfold\_search \
 --threads $SLURM\_CPUS\_PER\_TASK \
 for\_colabfold.fa $COLABFOLD\_DB pulldown**
[user@cn3144]$ **cd pulldown**
[user@cn3144]$ **ls -lh**
total 86M
-rw-r--r-- 1 user group  7.1M Sep 12 16:29 0.a3m
-rw-r--r-- 1 user group  9.2M Sep 12 16:29 10.a3m
-rw-r--r-- 1 user group  803K Sep 12 16:29 11.a3m
-rw-r--r-- 1 user group  611K Sep 12 16:29 12.a3m
[...snip...]
[user@cn3144]$ ## rename the alignments to something more useful
[user@cn3144]$ **rename\_colab\_search\_a3m.py**
Renaming 18.a3m to P62945.a3m
Renaming 20.a3m to Q9NX20.a3m
Renaming 0.a3m to P78344.a3m
Renaming 1.a3m to Q14240.a3m
[...snip...]
[user@cn3144]$ **ls -lh**
total 86M
-rw-r--r-- 1 user group  1.8M Sep 12 16:29 O76094.a3m
-rw-r--r-- 1 user group  9.2M Sep 12 16:29 P08240.a3m
-rw-r--r-- 1 user group  803K Sep 12 16:29 P09132.a3m
-rw-r--r-- 1 user group  2.1M Sep 12 16:29 P22090.a3m
[...snip...]
[user@cn3144]$ **cd ..**
[user@cn3144]$ **create\_individual\_features.py \
 --fasta\_paths=bait.fasta,candidates.fasta \
 --output\_dir=pulldown \
 --use\_precomputed\_msas=True \
 --max\_template\_date=2023-01-01 \
 --use\_mmseqs2=True \
 --skip\_existing=False**
[user@cn3144]$ ## exiting the sinteractive allocation using CPUs
[user@cn3144]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf]$

```

The local colabfold alignment pipleine above finished in 2.5h.


Notes on the alignment step:


* we do not recommend using `--use_mmseqs2=True` without pre-computing alignments
 with the local colabfold pileine.
* as mentioned in the [colabfold docs](https://hpc.nih.gov/apps/colabfold.html)
`colabfold_search` works more efficiently the more sequences are provided in a single
 query. Therefore all sequences for a pulldown should be submitted as a single job.


### Pairwise inference of multimer structures


This stage of the process will require a GPU. On a A100 GPU the full example
with 2 predictions per model woud take ≈ 17 h (1 bait \* 20 candidates \*
5 models/interaction \* 2 predictions/model \* 300s/prediction). Therefore we
only run the first 2 candidates in the example below. The output for this step
in this example is about 8.6GB/pair (172GB for all 20 pairs). If we had run
only 1 prediction per model per pair it would be half that.



```

[user@biowulf]$ **sinteractive --mem=50g -c8 --gres=lscratch:50,gpu:a100:1**
[user@cn3144]$ **module load alphapulldown**
[user@cn3144]$ **head -2 candidates.txt > 2candidates.txt** 
[user@cn3144]$ ## using more predictions / model will proportionally increase runtime
[user@cn3144]$ **run\_multimer\_jobs.py \
 --mode=pulldown \
 --num\_cycle=3 \
 --num\_predictions\_per\_model=2 \
 --output\_path=pulldown\_models \
 --protein\_lists=bait.txt,2candidates.txt \
 --monomer\_objects\_dir=pulldown**
[user@cn3144]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf]$

```

### Evaluation and visualization


alphapulldown has a script that will collect metrics from a run and create a jupyter notebook. Starting
this interactive session with `--tunnel` so that we can immediately look at the
notebook. See our [Jupyter documentation](https://hpc.nih.gov/apps/jupyter.html)
on how to connect to the notebook.



```

[user@biowulf]$ **sinteractive --tunnel --cpus-per-task=2 --mem=10g**
[user@cn3144]$ **module load alphapulldown**
[user@cn3144]$ **create\_notebook.py --help**
       USAGE: /opt/conda/envs/app/bin/create_notebook.py [flags]
flags:

/opt/conda/envs/app/bin/create_notebook.py:
  --[no]create_notebook: Whether creating a notebook
    (default: 'true')
  --cutoff: cutoff value of PAE. i.e. only pae<cutoff is counted good
    (default: '5.0')
    (a number)
  --output_dir: directory where predicted models are stored
    (default: '.')
  --pae_figsize: figsize of pae_plot, default is 50
    (default: '50')
    (an integer)
  --surface_thres: surface threshold. must be integer
    (default: '2')
    (an integer)

[user@cn3144]$ **cd pulldown\_models**
[user@cn3144]$ **create\_notebook.py**
[user@cn3144]$ **alphapulldown\_jupyter notebook --port $PORT1 --ip localhost --no-browser**

```

See [this](alphapulldown_notebook.html) partially rendered example notebook.


In addition alphapulldown provides a script to create a tabular summary
of putative hits:



```

[user@cn3144]$ **run\_get\_good\_pae.sh --help**

       USAGE: /app/alpha-analysis/get_good_inter_pae.py [flags]
flags:

/app/alpha-analysis/get_good_inter_pae.py:
  --cutoff: cutoff value of PAE. i.e. only pae<cutoff is counted good
    (default: '5.0')
    (a number)
  --output_dir: directory where predicted models are stored
  --surface_thres: surface threshold. must be integer
    (default: '2')
    (an integer)

[user@cn3144]$ **run\_get\_good\_pae.sh --output\_dir .**
[...much output...]
[user@cn3144]$ **head -2 predictions\_with\_good\_interpae.csv**
jobs,interface,Num_intf_residues,Polar,Hydrophobhic,Charged,contact_pairs, sc, hb, sb, int_solv_en, int_area,pi_score,pdb, pvalue,iptm_ptm,iptm,mpDockQ/pDockQ
P78344_and_P60842,C_B,100,0.23,0.25,0.39,104,0.52,32,27,-8.62,3795.28,0.88,,,0.7768986939967186,0.8371253,0.6383528254038463
[user@cn3144]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf]$

```

pDockQ and Protein Interface score (pi\_score) from the output table above



![](/images/alphapulldown.png)

### Custom mode inference


Custom mode can be used to explicitly define which interactions to test and can be used to only analyze fragments
of a sequence. This can be used, for example, to screen for interaction domains in proteins within very large proteins.
In our example we will use one of the highly scoring pairs above and map the interaction between the full length bait
and three fragments of Q14240 to illustrate the procedure. Running from the same directory as the examples above.



```

[user@biowulf]$ **sinteractive --mem=50g -c8 --gres=lscratch:50,gpu:a100:1**
[user@cn3144]$ **module load alphapulldown**
[user@cn3144]$ **cat <<\_\_EOF\_\_ > custom.txt
P78344;Q14240,1-130
P78344;Q14240,131-260
P78344;Q14240,161-407
\_\_EOF\_\_
[user@cn3144]$ **run\_multimer\_jobs.py \
 --mode=custom \
 --num\_cycle=3 \
 --num\_predictions\_per\_model=2 \
 --output\_path=custom \
 --protein\_lists=custom.txt \
 --monomer\_objects\_dir=pulldown**
[user@cn3144]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf]$**
```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
The alphafold based alignment and the structural inference steps can be parallelized using
slurm job arrays. The colabfold based alignment should submit *all* proteins in a single
alignment job.
### Alignment


#### Colabfold based pipeline - single job



```

#! /bin/bash

set -e

module load colabfold alphapulldown

cp ${ALPHAPULLDOWN_TEST_DATA:-none}/* .
cat bait.fasta candidates.fasta > for_colabfold.fa 

colabfold_search \
    --threads $SLURM_CPUS_PER_TASK \
    for_colabfold.fa $COLABFOLD_DB pulldown
pushd pulldown
rename_colab_search_a3m.py
popd

create_individual_features.py \
    --fasta_paths=bait.fasta,candidates.fasta \
    --output_dir=pulldown \
    --use_precomputed_msas=True \
    --max_template_date=2023-01-01 \
    --use_mmseqs2=True \
    --skip_existing=False

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```

[user@biowulf]$ **sbatch --cpus-per-task=16 --mem=128g --time=6:00:00 \
 --gres=lscratch:50 alphapulldown\_1a.sh**

```

#### Alphafold based pipeline - at most 10 concurrent jobs



```

#! /bin/bash

set -e

module load alphapulldown

cp ${ALPHAPULLDOWN_TEST_DATA:-none}/* .
create_individual_features.py \
    --fasta_paths=bait.fasta,candidates.fasta \
    --save_msa_files=True \
    --output_dir=pulldown \
    --use_precomputed_msas=False \
    --max_template_date=2023-01-01 \
    --skip_existing=True \
    --seq_index=$SLURM_ARRAY_TASK_ID

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```

[user@biowulf]$ **nbaits=$(grep -c ">" baits.fasta)**
[user@biowulf]$ **ncandidates=$(grep -c ">" candidates.fasta)**
[user@biowulf]$ **ntotal=$(( $nbaits + $ncandidates ))**
[user@biowulf]$ # Run a job array, 10 jobs at a time:
[user@biowulf]$ **sbatch --cpus-per-task=8 --mem=100g --time=2:00:00 --gres=lscratch:50 \
 --array=1-$ntotal%10 alphapulldown\_1b.sh**

```

### Pairwise inference of multimer structures - at most 10 concurrent jobs



```

#! /bin/bash

set -e

module load alphapulldown

run_multimer_jobs.py \
    --mode=pulldown \
    --num_cycle=3 \
    --num_predictions_per_model=2 \
    --output_path=pulldown_models \
    --protein_lists=bait.txt,candidates.txt \
    --monomer_objects_dir=pulldown \
    --job_index=$SLURM_ARRAY_TASK_ID

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```

[user@biowulf]$ **nbaits=$(wc -l < bait.txt)**
[user@biowulf]$ **ncandidates=$(wc -l < candidates.txt)**
[user@biowulf]$ **ntotal=$(( $nbaits \* $ncandidates ))**
[user@biowulf]$ # Run a job array, 10 jobs at a time:
[user@biowulf]$ **sbatch --cpus-per-task=8 --mem=50g --time=3:00:00 \
 --gres=lscratch:50,gpu:a100:1 --partition=gpu --array=1-$ntotal%10 alphapulldown\_2.sh**

```











