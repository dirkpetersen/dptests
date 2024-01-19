

document.querySelector('title').textContent = 'alphafold2 on Biowulf';

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

alphafold2 on Biowulf


|  |
| --- |
| 
Quick Links
[Changelog](#changes)
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int)
[Benchmarking](#bench)
[Batch job](#sbatch) 
[ColabFold](#colabfold)
 |




From the official documentation




> This package provides an implementation of the inference pipeline
>  of AlphaFold v2.0. This is a completely new model that was entered in
>  CASP14 and published in Nature. For simplicity, we refer to this model as
>  AlphaFold throughout the rest of this document.


### References:


* Jumper et al. *Highly accurate protein structure prediction with AlphaFold.*
 Nature 2021. [PubMed](https://pubmed.ncbi.nlm.nih.gov/34265844/)  | 
 [Journal](https://www.nature.com/articles/s41586-021-03819-2)



Changelog
[top](#top)

2023-05-23: alphafold 2.3.2 becomes the default version.
* No changes to model parameters or databases.
* `--run_relax` has been replaced with `--models_to_relax` with a default
 of `best`. That means that only the best model (ranked\_0.pdb) will be relaxed.
* The `jax` python module is no loger required to read the pickled results files

2023-05-10: updated mgnify to 2023-02
2023-02-13: alphafold 2.3.1 becomes the default version.
2023-02-10: alphafold 2.3.1 available.
Some notable changes (See [release page](https://github.com/deepmind/alphafold/releases) for full details):* Added new AlphaFold-Multimer models with better accuracy on large protein complexes.
* Added early stopping to recycling.
* Used bf16 in multimer inference – reduces GPU memory usage.
* Relaxation metrics are now saved in relax\_metrics.json.

2023-02-09: In place database update. This will apply to all alphafold version with the exception
 of the parameters which are version specific.
* uniref90 - 2022\_05
* mgnify - 2022\_05
* uniref30 - 2021\_06 (unchanged since last update; this used to be called uniclast30)
* bfd - casp14 (unchanged since last update)
* uniprot - 2022\_05
* pdb70 - 220313
* pdb mmcif - 230208
* pdb\_seqres - 230208



2022-09-21: the `--hhblits_extra_opts` option was ported from msa to run\_singularity
In a small number of cases hhblits fails to create alignments. This option can be used
 to fine tune the hhblits run (see below). Example: `run_singularity --hhblits_extra_opts="-maxres 80000 -prepre_smax_thresh 50" ...`
2022-07-11: the `msa` utility script has been disabled
Large scale use of the `msa` script may have been implicated in file system problems.
 The script has been removed until futher notice.
2022-06-02: added alphapickle to alphafold 2.2.0.
[alphapickle](https://github.com/mattarnoldbio/alphapickle) will be included in alphafold
 installs ≥ 2.2.0
2022-04-22: Version 2.2.0 becomes the default
* See [official announcement](https://github.com/deepmind/alphafold/releases/tag/v2.2.0)
* Switched to new AlphaFold-Multimer models with reduced numbers of clashes on average
 and slightly increased accuracy.
* removed `--is_prokaryote_list` option from `run_singularity`
 and `msa` as the prokaryotic pairing algorithm did not actually improve
 accuracy on average
* Added `--num_multimer_predictions_per_model=N` option to
 `run_singularity`. Runs N predictions per *multimer* model - each
 with different seeds. Defaults to 5 and will increase runtime
* Added `--model_config` option to `run_singularity`. This
 allows users to use a customized `alphafold/model/config.py` to alphafold
 for tuning certain parameters



2022-02-22: Version 2.1.2 becomes the default
* No changes to network weights
* Relaxation now defaults to using GPU
* Added options for controlling relaxation: `--[no]run_relax` and `--[no]enable_gpu_relax`
* New script to generate the multiple sequence alignments (MSAs) only (`msa`). This is 
 highly recommended as MSA generation consumes 
 about 60% of runtime and does not make use of GPU. The workflow is to run msa generation on CPU nodes
 and then run model prediction on a GPU node with `run_singularity --use_precomputed_msas ...`.
 This script also supports adding extra options to hhblits or overriding some alphafold defaults with option
 `--hhblits_extra_opts`.
 Note that this script cannot run on x2650 nodes because it depends on an AVX2 hhblits build.

2021-11-15: Version 2.1.1 becomes the default
* the `--preset` option was renamed to --model\_preset
* the `--db_preset=reduced_dbs` option is now supported for all alphafold versions
* The inofficial `--use_ptm` option became obsolete with introduction of the official `--model_preset`
 option and was removed
* This version uses the new model parameters (2021-10-27) released with 2.1.1. This includes
 the parameters for the multimer model. The 2.0.X modules will continue to use the previous parameters.

2021-11-14: Database update
Databases were updated in place: pdb mmcif and pdb70 (211110). New databases only used by multimer model: pdb\_seqres,
 uniprot
2021-10-19: Added `--use_ptm` option to run\_singularity
Use the pTM models, which were fine-tuned to produce pTM (predicted TM-score) and predicted 
 aligned error values alongside their structure predictions.
2021-10-18: Adaptation of the alphafold\_advanced notebook from [ColabFold](https://github.com/sokrypton/ColabFold)
 available in version 2.0.1.
Allows prediction of protein complexes with unmodified alphafold network weights. So far only an interactive notebook is available. See
 [below](#colabfold) for more details
2021-10-01: Version 2.0.0-24-g1d43aaf was tagged as 2.0.1
The modules for 2.0.0-24-g1d43aaf and 2.0.1 point to the same installation since
 the release was tagged after this revision was installed.
2021-09-21: Version 2.0.0-24-g1d43aaf becomes the default version on biowulf
Most noticable change should be the inclusion of pLDDT in the PDB B-factor column
2021-09-16: Database update (in place)
The following databases used by alphafold were updated in place: mgnify (2018\_12 to 2019\_05),
 pdb70 (200401 to 210901), pdb mmcif (210717 to 210915, 1969 additional structures), uniclust30 (2018\_08 to
 2021\_06 from <http://gwdu111.gwdg.de/~compbiol/uniclust/2021_06/>).
 Uniref90 and BFD are unchanged.

Documentation
[top](#top)
* alphafold2 on [GitHub](https://github.com/deepmind/alphafold)
* [Blog post](https://deepmind.com/blog/article/alphafold-a-solution-to-a-50-year-old-grand-challenge-in-biology)


Important Notes
* Module Name: alphafold2 (see [the modules page](/apps/modules.html) for more information)
* Alphafold2 first runs some multithreaded analyses using up to 8 CPUs before running
 model inference on the GPU. At this point these steps can't be separated and therefore
 for the first step of the job the GPU will remain idle.
* Example files in `$ALPHAFOLD2_TEST_DATA`
* Reference data in `/fdb/alphafold2/`
* Alphafold2 expects input to be *upper case* amino acid sequences



Interactive job

[top](#top)
Allocate an [interactive session](/docs/userguide.html#int) and run the program. In this example the whole pipeline
including multiple sequence alignment and model predictions are run with `run_singularity` on a GPU node.



```

[user@biowulf]$ **sinteractive --mem=60g --cpus-per-task=8 --gres=lscratch:100,gpu:v100x:1**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144]$ **module load alphafold2/2.3.1**

```

To predict the structure of a protein already in PDB without using its
experimental structure as a template set `max_template_date` to
before the release date of the structure. For example, to reproduce the T1049
CASP14 target with 144 aa. On a V100x this prediction runs for about 1h.



```

[user@cn3144]$ **run\_singularity --helpfull**  # use --help for shorter help message
Singularity launch script for Alphafold.
flags:
/usr/local/apps/alphafold2/2.3.1/bin/run_singularity:
  --[no]benchmark: Run multiple JAX model evaluations to obtain a timing that excludes the compilation time,
    which should be more indicative of the time required for inferencing many proteins.
    (default: 'false')
  --db_preset: <full_dbs|reduced_dbs>: Choose preset MSA database configuration - smaller genetic
    database config (reduced_dbs) or full genetic database config (full_dbs)
    (default: 'full_dbs')
  --[no]dry_run: Print command that would have been executed and exit.
    (default: 'false')
  --[no]enable_gpu_relax: Run relax on GPU if GPU is enabled.
    (default: 'true')
  --fasta_paths: Paths to FASTA files, each containing a prediction target that will be folded one after
    another. If a FASTA file contains multiple sequences, then it will be folded as a multimer. Paths should
    be separated by commas. All FASTA paths must have a unique basename as the basename is used to name the
    output directories for each prediction. (a comma separated list)
  --gpu_devices: Comma separated list of devices to pass to NVIDIA_VISIBLE_DEVICES.
    (default: 'all')
  --max_template_date: Maximum template release date to consider (ISO-8601 format: YYYY-MM-DD). Important
    if folding historical test sets.
  --model_config: Use this file instead of default alphafold/model/config.py
  --model_preset: <monomer|monomer_casp14|monomer_ptm|multimer>: Choose preset model configuration -
    the monomer model, the monomer model with extra ensembling, monomer model with pTM head, or multimer model
    (default: 'monomer')
  --num_multimer_predictions_per_model: How many predictions (each with a different random seed) will be
    generated per model. E.g. if this is 2 and there are 5 models then there will be 10 predictions per
    input. Note: this FLAG only applies if model_preset=multimer
    (default: '5')
    (an integer)
  --output_dir: Path to a directory that will store the results.
  --[no]run_relax: Whether to run the final relaxation step on the predicted models. Turning relax off might
    result in predictions with distracting stereochemical violations but might help in case you are having
    issues with the relaxation stage.
    (default: 'true')
  --[no]use_gpu: Enable NVIDIA runtime to run with GPUs.
    (default: 'true')
  --[no]use_precomputed_msas: Whether to read MSAs that have been written to disk instead of running the
    MSA tools. The MSA files are looked up in the output directory, so it must stay the same between multiple
    runs that are to reuse the MSAs. WARNING: This will not check if the sequence, database or configuration
    have changed.
    (default: 'false')
...
absl.logging:
  --[no]alsologtostderr: also log to stderr?
    (default: 'false')
  --log_dir: directory to write logfiles into
    (default: '')
  --logger_levels: Specify log level of loggers. The format is a CSV list of `name:level`. Where `name` is the
    logger name used with `logging.getLogger()`, and `level` is a level name  (INFO, DEBUG, etc). e.g.
    `myapp.foo:INFO,other.logger:DEBUG`
    (default: '')
  --[no]logtostderr: Should only log to stderr?
    (default: 'false')
  --[no]showprefixforinfo: If False, do not prepend prefix to info messages when it's logged to stderr, --verbosity
    is set to INFO level, and python logging is used.
    (default: 'true')
  --stderrthreshold: log messages at this level, or more severe, to stderr in addition to the logfile.
    Possible values are 'debug', 'info', 'warning', 'error', and 'fatal'.  Obsoletes --alsologtostderr.
    Using --alsologtostderr cancels the effect of this flag. Please also note that this flag is
    subject to --verbosity and requires logfile not be stderr.
    (default: 'fatal')
  -v,--verbosity: Logging verbosity level. Messages logged at this level or lower will be included. Set to 1
    for debug logging. If the flag was not set or supplied, the value will be changed from the default of
    -1 (warning) to 0 (info) after flags are parsed.
    (default: '-1')
    (an integer)
...


[user@cn3144]$ **run\_singularity \
 --model\_preset=monomer \
 --fasta\_paths=$ALPHAFOLD2\_TEST\_DATA/T1049.fasta \
 --max\_template\_date=2022-12-31 \
 --output\_dir=$PWD**
[user@cn3144]$ **tree T1049**
T1049/
├── [user   1.1M]  features.pkl
├── [user   4.0K]  msas
│   ├── [user    33K]  bfd_uniclust_hits.a3m
│   ├── [user    18K]  mgnify_hits.sto
│   └── [user   121K]  uniref90_hits.sto
├── [user   170K]  ranked_0.pdb               # <-- shown below
├── [user   170K]  ranked_1.pdb
├── [user   170K]  ranked_2.pdb
├── [user   171K]  ranked_3.pdb
├── [user   170K]  ranked_4.pdb
├── [user    330]  ranking_debug.json
├── [user   170K]  relaxed_model_1.pdb
├── [user   170K]  relaxed_model_2.pdb
├── [user   170K]  relaxed_model_3.pdb
├── [user   170K]  relaxed_model_4.pdb
├── [user   171K]  relaxed_model_5.pdb
├── [user    11M]  result_model_1.pkl
├── [user    11M]  result_model_2.pkl
├── [user    11M]  result_model_3.pkl
├── [user    11M]  result_model_4.pkl
├── [user    11M]  result_model_5.pkl
├── [user    771]  timings.json
├── [user    87K]  unrelaxed_model_1.pdb
├── [user    87K]  unrelaxed_model_2.pdb
├── [user    87K]  unrelaxed_model_3.pdb
├── [user    87K]  unrelaxed_model_4.pdb
└── [user    87K]  unrelaxed_model_5.pdb


```

The processes prior to model inference on the GPU consumed up to 40 GB of
memory for this protein. Memory requirements will vary with different size
proteins.




![alignment for T1049 predicted and experimental structures](/images/alphafold_T1049.gif)
**Figure 1**. This is the highest confidence
 prediction (`ranked_0.pdb`, blue) aligned with the actual
 structure for this protein
 ([6Y4F](https://www.rcsb.org/structure/6y4f), green)


Note that the model .pkl files which, unlike the .pdb files, are not
re-ordered into ranked\_ files contain a lot of information about the models. These are
python pickle files and python can be used to explore and visualize them. For example:



```

[user@cn3144]$ **conda activate my\_py39** # needs jupyter and the packages imported below
[user@cn3144]$ **cd T1049**
[user@cn3144]$ **jupyter console**
In [1]: import pickle
In [2]: import json
In [3]: import pprint
In [4]: import jax   # only needed for version 2.3.1
In [5]: pprint.pprint(json.load(open("ranking_debug.json", encoding="ascii")))
{'order': ['model_2_pred_0',
           'model_3_pred_0',
           'model_1_pred_0',
           'model_5_pred_0',
           'model_4_pred_0'],
 'plddts': {'model_1_pred_0': 88.44386138278787,
            'model_2_pred_0': 91.83564104655056,
            'model_3_pred_0': 88.49961929441032,
            'model_4_pred_0': 86.73066329994059,
            'model_5_pred_0': 87.4009420322368}}
### so model 2 is the best model in this run and corresponds to ranked_0.pdf
In [6]: best_model = pickle.load(open("result_model_2_pred_0.pkl", "rb"))
In [7]: list(best_model.keys())
Out[7]:
['distogram',
 'experimentally_resolved',
 'masked_msa',
 'predicted_lddt',
 'structure_module',
 'plddt',
 'ranking_confidence']
In [8]: best_model['plddt'].shape
Out[8]: (141,)

```

The predicted alignment error (PAE) is only produced by the monomer\_ptm and
multimer models. Since version 2.2.0 we also include [alphapickle](https://github.com/mattarnoldbio/alphapickle)
with alphafold to create plots, csv files, and chimera attribute files for each ranked model. By default output
will be saved to the same folder. See `-h` for more options.



```

[user@cn3144]$ **alphapickle -od T1049**

```

If the model above was created with the monomer\_ptm model the following two plots are generated for each model:



![](/images/alphafold_t1049_ptm_alphapickle.png)
The next example shows how to run a multimer model (available from version
2.1.1). The example used is a recently published [PI3K structure](https://www.rcsb.org/structure/7MYN).



```

[user@cn3144]$ **cat $ALPHAFOLD2\_TEST\_DATA/pi3k.fa**
>sp|P27986|P85A_HUMAN Phosphatidylinositol 3-kinase regulatory subunit alpha OS=Homo sapiens OX=9606 GN=PIK3R1 PE=1 SV=2
MSAEGYQYRALYDYKKEREEDIDLHLGDILTVNKGSLVALGFSDGQEARPEEIGWLNGYN
ETTGERGDFPGTYVEYIGRKKISPPTPKPRPPRPLPVAPGSSKTEADVEQQALTLPDLAE
QFAPPDIAPPLLIKLVEAIEKKGLECSTLYRTQSSSNLAELRQLLDCDTPSVDLEMIDVH
VLADAFKRYLLDLPNPVIPAAVYSEMISLAPEVQSSEEYIQLLKKLIRSPSIPHQYWLTL
QYLLKHFFKLSQTSSKNLLNARVLSEIFSPMLFRFSAASSDNTENLIKVIEILISTEWNE
RQPAPALPPKPPKPTTVANNGMNNNMSLQDAEWYWGDISREEVNEKLRDTADGTFLVRDA
STKMHGDYTLTLRKGGNNKLIKIFHRDGKYGFSDPLTFSSVVELINHYRNESLAQYNPKL
DVKLLYPVSKYQQDQVVKEDNIEAVGKKLHEYNTQFQEKSREYDRLYEEYTRTSQEIQMK
RTAIEAFNETIKIFEEQCQTQERYSKEYIEKFKREGNEKEIQRIMHNYDKLKSRISEIID
SRRRLEEDLKKQAAEYREIDKRMNSIKPDLIQLRKTRDQYLMWLTQKGVRQKKLNEWLGN
ENTEDQYSLVEDDEDLPHHDEKTWNVGSSNRNKAENLLRGKRDGTFLVRESSKQGCYACS
VVVDGEVKHCVINKTATGYGFAEPYNLYSSLKELVLHYQHTSLVQHNDSLNVTLAYPVYA
QQRR
>sp|P42336|PK3CA_HUMAN Phosphatidylinositol 4,5-bisphosphate 3-kinase catalytic subunit alpha isoform OS=Homo sapiens OX=9606 GN=PIK3CA PE=1 SV=2
MPPRPSSGELWGIHLMPPRILVECLLPNGMIVTLECLREATLITIKHELFKEARKYPLHQ
LLQDESSYIFVSVTQEAEREEFFDETRRLCDLRLFQPFLKVIEPVGNREEKILNREIGFA
IGMPVCEFDMVKDPEVQDFRRNILNVCKEAVDLRDLNSPHSRAMYVYPPNVESSPELPKH
IYNKLDKGQIIVVIWVIVSPNNDKQKYTLKINHDCVPEQVIAEAIRKKTRSMLLSSEQLK
LCVLEYQGKYILKVCGCDEYFLEKYPLSQYKYIRSCIMLGRMPNLMLMAKESLYSQLPMD
CFTMPSYSRRISTATPYMNGETSTKSLWVINSALRIKILCATYVNVNIRDIDKIYVRTGI
YHGGEPLCDNVNTQRVPCSNPRWNEWLNYDIYIPDLPRAARLCLSICSVKGRKGAKEEHC
PLAWGNINLFDYTDTLVSGKMALNLWPVPHGLEDLLNPIGVTGSNPNKETPCLELEFDWF
SSVVKFPDMSVIEEHANWSVSREAGFSYSHAGLSNRLARDNELRENDKEQLKAISTRDPL
SEITEQEKDFLWSHRHYCVTIPEILPKLLLSVKWNSRDEVAQMYCLVKDWPPIKPEQAME
LLDCNYPDPMVRGFAVRCLEKYLTDDKLSQYLIQLVQVLKYEQYLDNLLVRFLLKKALTN
QRIGHFFFWHLKSEMHNKTVSQRFGLLLESYCRACGMYLKHLNRQVEAMEKLINLTDILK
QEKKDETQKVQMKFLVEQMRRPDFMDALQGFLSPLNPAHQLGNLRLEECRIMSSAKRPLW
LNWENPDIMSELLFQNNEIIFKNGDDLRQDMLTLQIIRIMENIWQNQGLDLRMLPYGCLS
IGDCVGLIEVVRNSHTIMQIQCKGGLKGALQFNSHTLHQWLKDKNKGEIYDAAIDLFTRS
CAGYCVATFILGIGDRHNSNIMVKDDGQLFHIDFGHFLDHKKKKFGYKRERVPFVLTQDF
LIVISKGAQECTKTREFERFQEMCYKAYLAIRQHANLFINLFSMMLGSGMPELQSFDDIA
YIRKTLALDKTEQEALEYFMKQMNDAHHGGWTTKMDWIFHTIKQHALN

[user@cn3144]$ **run\_singularity \
 --fasta\_paths=$ALPHAFOLD2\_TEST\_DATA/pi3k.fa \
 --max\_template\_date=2021-11-01 \
 --model\_preset multimer \
 --num\_multimer\_predictions\_per\_model=2 \
 --output\_dir=$PWD**
...snip...
[user@cn3144]$ **exit**

```



![PI3K heterodimer model animation](/images/alphafold_pi3k.gif)
**Figure 2**. Best alphafold
 model for Phosphoinositide 3-kinase alpha (PI3Kα) model obtained
 in the example above. The two subunits are shown in blue (catalytic
 subunit, p110) and green (regulatory subunit, p85), respectively and
 shaded by pLDDT from light (low) to dark (high). Comparision with the
 Cryo-EM structure ([7MYN](https://www.rcsb.org/structure/7MYN))
 showed close agreement and some high confidence
 predicitons for areas that did not resolve in the
 published structure. 


Benchmarking
[top](#top)
To get an idea of runtimes of alphafold2 we first ran 4 individual proteins
on all our available GPUs. The proteins ranged in size from 144 aa to 622 aa. Note
that for all but the smallest protein, K80 GPUs were not suitable and should not
be considered for alphafold2. These tests were run with default settings except for
a fixed `--max_template_date=2021-07-31`




![alphafold2 benchmark plots](/images/alphafold2_benchmarks.png)
**Figure 3**. Alphafold runtimes for 4 proteins with 8 CPUs and 60GB of
 memory. T1049: CASP14 target,
 144 aa. D10R: Vaccinia virus WR protein D10R, 248 aa. A11R: Vaccinia virus WR protein A11R,
 318 aa. T1036s1: CASP14 target, 622 aa. Note that (1) k80 GPUs are not suitable (2)
 2 GPUs for single proteins don't reduce runtime (3) p100s are about as good as the more modern
 GPUs for these examples (4) Runtimes were quite variable. (5) The degree
 to which the 8 CPUs were overloaded depended on protein size and overloading
 appeared to be most severe during the relaxation phase.


The runtime to run all 4 protein on a V100x GPU with 8 CPUs and 60GB of memory was 3.2h, 
slightly less than the individual runtimes of the 4 proteins run separately. For this one job we also
increased the number of CPUs to 16 or the number of GPUs to 2, neither of which appeared to shorted
the runtime


The resource usage profile of the combined alphafold2 pipeline in our testing thus far
is suboptimal and variable. Steps probably should be segregated into individual jobs with
proper resources. We hope to optimize this in the future


Batch job

[top](#top)
Note: when running multiple alphafold predictions please use the msa script available for
alphafold >=2.1.2 to precompute the multiple sequence alignments on CPU nodes and use GPU only for
model predictions. This is shown in the example below as 2 batch jobs.



```

#!/bin/bash
module load alphafold2/2.2.0
run_singularity \
    --model_preset=monomer \
    --fasta_paths=$ALPHAFOLD2_TEST_DATA/T1049.fasta \
    --max_template_date=2020-05-14 \
    --output_dir=$PWD

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=6 --partition=gpu --mem=40g --gres=gpu:v100x:1,lscratch:100 alphafold2_model.sh
```

ColabFold alphafold2\_advanced notebook

[top](#top)
This is an adaptation of the Alphafold2\_advanced
[Colabfold](https://github.com/sokrypton/ColabFold) notebook to
biowulf. Currently only the notebook is available but a batch wrapper is
in development.


Notes


* Currently only available for version alphafold 2.0.1
* jackhmmer is not yet implemented. The notebook currently only
 uses mmseqs2 via the online API to create the input MAS
* Like the colabfold notebook this does *not* use templates
 and uses unmodified alphafold weights trained on single proteins.
* Colab specific input forms had to be removed. Instead, input has
 to be provided in normal cells with options described in the text


The jupyter function will use an existing tunnel if it has been set up
with the `--tunnel` option to `sinteractive`. If there is no pre-existing tunnel, it will attempt to
set one up itself. That means it is possible to start the jupyter server in a
batch job and obtain the command to set up the tunnel from your computer to
the login node from the batch output. See our [tunneling documentation](/docs/tunneling/) for more information


Example use



```

[user@biowulf ~]$ **sinteractive --gres=lscratch:20,gpu:v100x:1 -c16 --mem=60g --tunnel**
salloc.exe: Pending job allocation 25316671
salloc.exe: job 25316671 queued and waiting for resources
salloc.exe: job 25316671 has been allocated resources
salloc.exe: Granted job allocation 25316671
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3299 are ready for job
srun: error: x11: no local DISPLAY defined, skipping
error: unable to open file /tmp/slurm-spank-x11.25316671.0
slurmstepd: error: x11: unable to read DISPLAY value

Created 1 generic SSH tunnel(s) from this compute node to
biowulf for your use at port numbers defined
in the $PORTn ($PORT1, ...) environment variables.


Please create a SSH tunnel from your workstation to these ports on biowulf.
On Linux/MacOS, open a terminal and run:

    ssh  -L 46243:localhost:46243 user@biowulf.nih.gov

    For Windows instructions, see https://hpc.nih.gov/docs/tunneling

[user@cn3144]$ **module load alphafold2/2.0.1**
[user@cn3144]$ **af\_colabfold\_advanced jupyter name\_of\_notebook\_copy.ipynb**


```

The `af_colabfold_advanced` script creates a copy of the notebook
in your current directory that you can work with and starts the jupyter server
up. Once the provided ssh command is used to establish a tunnel to the login node
exactly as when using regular jupyter notebooks, the notebook can be visited
at the address provided by jupyter during startup






