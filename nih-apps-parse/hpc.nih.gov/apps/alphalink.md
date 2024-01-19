

document.querySelector('title').textContent = "alphalink";

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

AlphaLink on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
 |


 AlphaLink predicts protein structures using deep learning given a sequence and a set of experimental contacts. 
 It extends OpenFold with crosslinking MS data or other experimental distance restraint by explicitly incorporating
 them in the OpenFold architecture.


### References:


* "Protein structure prediction with in-cell photo-crosslinking mass spectrometry and deep learning", Nat. Biotech. XXX doi:10.1038/s41587-023-01704-z.
 [PubMed](https://pubmed.ncbi.nlm.nih.gov/36941363/) | 
 [Journal](https://www.nature.com/articles/s41587-023-01704-z)


Documentation
* AlphaLink on [GitHub](https://github.com/lhatsk/AlphaLink)


Important Notes
* Module Name: `alphalink` (see [the modules page](/apps/modules.html) 
 for more information)
* Some steps in an analysis require GPU
* Example files in `$ALPHALINK_TEST_DATA`
* Environment variables set:
	+ `ALPHALINK_CP_[CACA|DIST]` - AlphaLink model checkpoint files
	+ `[UNIREF90|MGNIFY|PDB70|UNICLUST30]_PATH` - Path to Alphafold DBs
	+ `[JACKHMMER|HHBLITS|HHSEARCH|KALIGN]_BIN` - Required programs



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
In the following example we will predict a structure from a single FASTA sequence and the
 crosslinking mass-spectrometry residue pairs.


### AlphaLink structure prediction


#### AlphaLink with CSV crosslink restraints



```

[user@biowulf]$ **sinteractive --gres=gpu:1,lscratch:100 --constraint='gpup100|gpuv100|gpuv100x|gpua100' -c 8 --mem=32g**
[user@cn3144]$ **module load alphalink/1.0**
[user@cn3144]$ **cd /lscratch/$SLURM\_JOB\_ID**
[user@cn3144]$ **cp -r ${ALPHALINK\_TEST\_DATA:-none}/\* .**
[user@cn3144]$ **time predict\_with\_crosslinks.py CDK/fasta/CDK.fasta CDK/crosslinks/mixture.csv \
 --checkpoint\_path $ALPHALINK\_CP\_CACA \
 --uniref90\_database\_path $UNIREF90\_PATH \
 --mgnify\_database\_path $MGNIFY\_PATH \
 --pdb70\_database\_path $PDB70\_PATH \
 --uniclust30\_database\_path $UNICLUST30\_PATH \
 --jackhmmer\_binary\_path $JACKHMMER\_BIN \
 --hhblits\_binary\_path $HHBLITS\_BIN \
 --hhsearch\_binary\_path $HHSEARCH\_BIN \
 --kalign\_binary\_path $KALIGN\_BIN**
INFO:/usr/local/apps/alphalink/1.0/AlphaLink/predict_with_crosslinks.py:Loaded OpenFold parameters at /fdb/alphalink/finetuning_model_5_ptm_CACA_10A.pt...
INFO:/usr/local/apps/alphalink/1.0/AlphaLink/predict_with_crosslinks.py:Generating alignments for sp|P24941|CDK2_HUMAN...
INFO:/usr/local/apps/alphalink/1.0/AlphaLink/predict_with_crosslinks.py:Loaded 9 restraints...
INFO:/usr/local/apps/alphalink/1.0/AlphaLink/predict_with_crosslinks.py:Running inference for sp|P24941|CDK2_HUMAN...
INFO:/usr/local/apps/alphalink/1.0/AlphaLink/predict_with_crosslinks.py:Inference time: 56.18929745070636
INFO:/usr/local/apps/alphalink/1.0/AlphaLink/predict_with_crosslinks.py:Output written to /lscratch/0000000/predictions/sp|P24941|CDK2_HUMAN_model_5_ptm_crosslinks_unrelaxed.pdb...
INFO:/usr/local/apps/alphalink/1.0/AlphaLink/predict_with_crosslinks.py:Relaxed output written to /lscratch/0000000/predictions/sp|P24941|CDK2_HUMAN_model_5_ptm_crosslinks_relaxed.pdb...
INFO:/usr/local/apps/alphalink/1.0/AlphaLink/predict_with_crosslinks.py:Model output written to /lscratch/0000000/predictions/sp|P24941|CDK2_HUMAN_model_5_ptm_crosslinks_output_dict.pkl...
real    44m12.596s
user    159m13.407s
sys     4m3.506s
[user@cn3144]$ **ls -lh predictions**
total 153M
-rw-r--r-- 1 user group 153M Oct 17 12:26 'sp|P24941|CDK2_HUMAN_model_5_ptm_crosslinks_output_dict.pkl'
-rw-r--r-- 1 user group 384K Oct 17 12:26 'sp|P24941|CDK2_HUMAN_model_5_ptm_crosslinks_relaxed.pdb'
-rw-r--r-- 1 user group 190K Oct 17 12:25 'sp|P24941|CDK2_HUMAN_model_5_ptm_crosslinks_unrelaxed.pdb'

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
### Prediction of structure with crosslinking constraints - at most 10 concurrent jobs



```

#! /bin/bash

set -e

module load alphalink

predict_with_crosslinks.py protein.fasta crosslinks.csv \
                    --checkpoint_path $ALPHALINK_CP_CACA \
                    --uniref90_database_path $UNIREF90_PATH \
                    --mgnify_database_path $MGNIFY_PATH \
                    --pdb70_database_path $PDB70_PATH \
                    --uniclust30_database_path $UNICLUST30_PATH \
                    --jackhmmer_binary_path $JACKHMMER_BIN \
                    --hhblits_binary_path $HHBLITS_BIN \
                    --hhsearch_binary_path $HHSEARCH_BIN \
                    --kalign_binary_path $KALIGN_BIN

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```

[user@biowulf]$ **sbatch --cpus-per-task=8 --mem=32g --time=3:00:00 \
 --gres=lscratch:50,gpu:1 \
 --partition=gpu \
 --constraint='gpup100|gpuv100|gpuv100x|gpua100' \
 alphalink.sh**

```







