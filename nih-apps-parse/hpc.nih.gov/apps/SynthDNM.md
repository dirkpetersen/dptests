

document.querySelector('title').textContent = 'SynthDNM: a random-forest based classifier for robust de novo prediction of SNPs and indels.';
**SynthDNM: a random-forest based classifier for robust de novo prediction of SNPs and indels.**


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



SynthDNM is a random-forest based classifier that can be readily adapted
to new sequencing or variant-calling pipelines 
by applying a flexible approach to constructing simulated
training examples from real data. The optimized SynthDNM classifiers predict de novo SNPs and indels
with robust accuracy across multiple methods of variant calling.



### References:


 Aojie Lian, James Guevara, Kun Xia and Jonathan Sebat   

 *Customized de novo mutation detection for any variant
calling pipeline: SynthDNM*    

[Bioinformatics](https://academic.oup.com/bioinformatics/article/37/20/3640/6209072), Volume 37, Issue 20, October 2021, Pages 3640–3641, https://doi.org/10.1093/bioinformatics/btab225


Documentation
* [SynthDNM Github page](https://github.com/james-guevara/synthdnm)


Important Notes
* Module Name: SynthDNM (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set
	+ **SDNM\_HOME**  installation directory
	+ **SDNM\_BIN**       executable directory
	+ **SDNM\_DATA**    sample data folder



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. 
  
Sample session on a GPU node:



```

[user@biowulf ~]$ **sinteractive --mem=8g -c4 --gres=lscratch:20**
[user@cn2379 ~]$ **module load synthdnm**
[+] Loading singularity  3.8.4  on cn2379
[+] Loading synthdnm  1.1.50

```

Basic usage:

```

[user@cn2379 ~]$ **run\_synthdnm.py -h**
usage: run_synthdnm.py [-h] [--vcf_file VCF_FILE] --ped_file PED_FILE
                       [--region REGION] [--features_file FEATURES_FILE]
                       [--output_folder OUTPUT_FOLDER]
                       [--training_set_tsv TRAINING_SET_TSV]
                       {classify,make_training_set,train,grid_search} ...

SynthDNM: a de novo mutation classifier and training paradigm

positional arguments:
  {classify,make_training_set,train,grid_search}
                        Available sub-commands
    classify            Classify DNMs using pre-trained classifiers.
    make_training_set   Make training set.
    train               Train classifiers
    grid_search         Randomized grid search across hyperparameters.

optional arguments:
  -h, --help            show this help message and exit
  --vcf_file VCF_FILE   VCF file input
  --ped_file PED_FILE   Pedigree file (.fam/.ped/.psam) input
  --region REGION       Interval ('{}' or '{}:{}-{}' in format of chr or
                        chr:start-end) on which to run training or
                        classification
  --features_file FEATURES_FILE
                        Features file input
  --output_folder OUTPUT_FOLDER
                        Output folder for output files (if not used, then
                        output folder is set to 'synthdnm_output')
  --training_set_tsv TRAINING_SET_TSV
                        Training set file (created using make_training_set
                        mode)
[user@cn2379 ~]$ **run\_synthdnm.py classify -h**
usage: run_synthdnm.py classify [-h] --clf_folder CLF_FOLDER
                                [-feature_extraction_only]

optional arguments:
  -h, --help            show this help message and exit
  --clf_folder CLF_FOLDER
                        Folder that contains the classifiers, which must be in
                        .pkl format (if not specified, will look for them in
                        the default data folder)
  -feature_extraction_only
                        Only output the features file (without classifying

```


End the interactive session:

```

[user@cn2379 ~]$ **exit**

```

Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. synthdnm.sh). For example:



```

#!/bin/bash
set -e
module load synthdnm
cp $SDNM_DATA/* .
synthdnm -v tutorial.vcf -f tutorial.ped

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
**sbatch synthdnm.sh**
```





