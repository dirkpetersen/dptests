

document.querySelector('title').textContent = 'TensorQTL: a GPU-enabled, ultrafast QTL mapper';
**TensorQTL: a GPU-enabled, ultrafast QTL mapper**


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



TensoorQTL leverages general-purpose libraries and graphics processing units
(GPUs) to achieve high efficiency of computations at low costR. Using PyTorch or TensorFlow
it allows > 200-fold decreases in runtime and ~ 5–10-fold reductions in cost when
running on GPUs relative to CPUs. 




### References:


* Amaro Taylor-Weiner, François Aguet, Nicholas J. Haradhvala, Sager Gosai, Shankara Anand, Jaegil Kim, Kristin Ardlie, Eliezer M. Van Allen and Gad Getz   

 *Scaling computational genomics to millions of individuals with GPUs.*   

[Genome Biology (2019) 20:228,](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1836-7)   
*


Documentation
* [TensorQTL Github page](https://github.com/broadinstitute/tensorqtl)


Important Notes
* Module Name: TensorQTL (see [the modules page](/apps/modules.html) for more information)
* Singlethreaded
* Unusual environment variables set
	+ **TQTL\_HOME**   installation directory
	+ **TQTL\_BIN**         executable directory
	+ **TQTL\_DATA**     test data directory* Example files in **$TQTL\_EXAMPLE**



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=24g --gres=gpu:p100:1,lscratch:50 -c8**
[user@cn4199 ~]$ **module load TensorQTL**
[+] Loading singularity  3.8.5-1 
[+] Loading cuDNN/7.6.5/CUDA-10.2 libraries...
[+] Loading CUDA Toolkit  10.2.89  ...
[+] Loading TensorQTL 1.0.7  ...

```

Usage: 

```

[user@cn4199 ~]$ **tensorqtl -h**
usage: tensorqtl [-h] [--mode {cis,cis_nominal,cis_independent,cis_susie,trans}] [--covariates COVARIATES]
                 [--permutations PERMUTATIONS] [--interaction INTERACTION] [--cis_output CIS_OUTPUT]
                 [--phenotype_groups PHENOTYPE_GROUPS] [--window WINDOW] [--pval_threshold PVAL_THRESHOLD]
                 [--maf_threshold MAF_THRESHOLD] [--maf_threshold_interaction MAF_THRESHOLD_INTERACTION] [--return_dense]
                 [--return_r2] [--best_only] [--output_text] [--batch_size BATCH_SIZE] [--load_split] [--warn_monomorphic]
                 [--fdr FDR] [--qvalue_lambda QVALUE_LAMBDA] [--seed SEED] [-o OUTPUT_DIR]
                 genotype_path phenotype_bed prefix

tensorQTL: GPU-based QTL mapper

positional arguments:
  genotype_path         Genotypes in PLINK format
  phenotype_bed         Phenotypes in BED format
  prefix                Prefix for output file names

optional arguments:
  -h, --help            show this help message and exit
  --mode {cis,cis_nominal,cis_independent,cis_susie,trans}
                        Mapping mode. Default: cis
  --covariates COVARIATES
                        Covariates file, tab-delimited, covariates x samples
  --permutations PERMUTATIONS
                        Number of permutations. Default: 10000
  --interaction INTERACTION
                        Interaction term(s)
  --cis_output CIS_OUTPUT
                        Output from 'cis' mode with q-values. Required for independent cis-QTL mapping.
  --phenotype_groups PHENOTYPE_GROUPS
                        Phenotype groups. Header-less TSV with two columns: phenotype_id, group_id
  --window WINDOW       Cis-window size, in bases. Default: 1000000.
  --pval_threshold PVAL_THRESHOLD
                        Output only significant phenotype-variant pairs with a p-value below threshold. Default: 1e-5 for
                        trans-QTL
  --maf_threshold MAF_THRESHOLD
                        Include only genotypes with minor allele frequency >= maf_threshold. Default: 0
  --maf_threshold_interaction MAF_THRESHOLD_INTERACTION
                        MAF threshold for interactions, applied to lower and upper half of samples
  --return_dense        Return dense output for trans-QTL.
  --return_r2           Return r2 (only for sparse trans-QTL output)
  --best_only           Only write lead association for each phenotype (interaction mode only)
  --output_text         Write output in txt.gz format instead of parquet (trans-QTL mode only)
  --batch_size BATCH_SIZE
                        Batch size. Reduce this if encountering OOM errors.
  --load_split          Load genotypes into memory separately for each chromosome.
  --warn_monomorphic    Warn if monomorphic variants are found.
  --fdr FDR             FDR for cis-QTLs
  --qvalue_lambda QVALUE_LAMBDA
                        lambda parameter for pi0est in qvalue.
  --seed SEED           Seed for permutations.
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Output directory

```

Running the test example:

```

[user@cn4199 ~]$ **cp $TQTL\_DATA/\* .**
[user@cn4199 ~]$ **cp $TQTL\_EXAMPLE/\*.py .**
[user@cn4199 ~]$ **python tensorqtl\_examples.py &**
python tensorqtl_examples.py &[1] 3000
PyTorch 1.12.0+cu102
Pandas 1.4.3

Mapping files:   0%|                                                                                     | 0/3 [00:00<?, ?it/s]
Mapping files: 100%|█████████████████████████████████████████████████████████████████████████████| 3/3 [00:23<00:00,  7.80s/it]
cis-QTL mapping: nominal associations for all variant-phenotype pairs
  * 445 samples
  * 301 phenotypes
  * 26 covariates

nvidia-smi
Mon Oct  3 10:14:00 2022
+-----------------------------------------------------------------------------+
| NVIDIA-SMI 470.82.01    Driver Version: 470.82.01    CUDA Version: 11.4     |
|-------------------------------+----------------------+----------------------+
| GPU  Name        Persistence-M| Bus-Id        Disp.A | Volatile Uncorr. ECC |
| Fan  Temp  Perf  Pwr:Usage/Cap|         Memory-Usage | GPU-Util  Compute M. |
|                               |                      |               MIG M. |
|===============================+======================+======================|
|   0  Tesla P100-PCIE...  On   | 00000000:0D:00.0 Off |                  Off |
| N/A   29C    P0    30W / 250W |    789MiB / 16280MiB |      1%      Default |
|                               |                      |                  N/A |
+-------------------------------+----------------------+----------------------+

+-----------------------------------------------------------------------------+
| Processes:                                                                  |
|  GPU   GI   CI        PID   Type   Process name                  GPU Memory |
|        ID   ID                                                   Usage      |
|=============================================================================|
|    0   N/A  N/A      3030      C   /usr/bin/python3                  787MiB |
+-----------------------------------------------------------------------------+
[user@cn2375 TensorQTL]$   * 13369268 variants
  * checking phenotypes: 301/301
  * Computing associations
    Mapping chromosome chr18
    processing phenotype 301/301
    time elapsed: 0.05 min
    * writing output
done.
cis-QTL mapping: empirical p-values for phenotypes
  * 445 samples
  * 301 phenotypes
  * 26 covariates
  * 13369268 variants
  * using seed 123456
  * checking phenotypes: 301/301
  * computing permutations
    processing phenotype 301/301
  Time elapsed: 0.40 min
done.
trans-QTL mapping
  * 445 samples
  * 19836 phenotypes
  * 26 covariates
  * 13369268 variants
    processing batch 1337/1337
    elapsed time: 0.62 min
  * 7620376 variants passed MAF >= 0.05 filtering
done.
[user@cn4199 ~]$ **exit**
salloc.exe: Relinquishing job allocation 59748321
[user@biowulf ~]$

```





