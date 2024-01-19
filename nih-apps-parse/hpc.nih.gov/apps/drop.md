

document.querySelector('title').textContent = 'drop on Biowulf';
drop on Biowulf


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



The drop is a pipeline to find aberrant events in RNA-Seq data, useful for diagnosis of rare disorders. 






### References:


* Yépez, V.A., Mertes, C., Müller, M.F. et al. *Detection of aberrant gene expression events in RNA sequencing data.*  Nat Protoc. 16, 1276–1296 (2021).
 [PubMed](https://pubmed.ncbi.nlm.nih.gov/33462443/) | 
 [Journal](https://www.nature.com/articles/s41596-020-00462-5)


Documentation
* drop Github:[Github](https://github.com/gagneurlab/drop)


Important Notes
* Module Name: drop (see [the modules page](/apps/modules.html) for more information)
 * Current drop command lines could be run as:
 
```

	drop
	
```
* The current version of drop was installed in it's own conda environment, so please do not load snakemake or R module.

Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive --cpus-per-task=10 --mem=10G**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load drop**
[user@cn3144 ~]$ **mkdir /data/$USER/drop\_test/**
[user@cn3144 ~]$ **cd /data/$USER/drop\_test/**
[user@cn3144 ~]$ **drop demo**
[user@cn3144 ~]$ **snakemake -n #dryrun**
[user@cn3144 ~]$ **snakemake --core 10**
WARNING: Using the mae defined genome instead of the globally defined one.
This will be deprecated in the future to allow for reference genomes to be defined in the sample annotation table. Please update your config and sample annotation table

WARNING: GENE_ANNOTATION must be a column in the sample annotation table, ANNOTATION is the old column name and will be deprecated in the future

WARNING: Less than 30 IDs in DROP_GROUP outrider
WARNING: Less than 30 IDs in DROP_GROUP import_exp
WARNING: Less than 30 IDs in DROP_GROUP fraser
check for missing R packages
Structuring dependencies...
Dependencies file generated at: /tmp/tmp1a00_tgp

Building DAG of jobs...
Using shell: /usr/bin/bash
Provided cores: 2
Rules claiming more threads will be scaled down.
Job counts:
	count	jobs
	1	AberrantExpression_Overview_R
	1	AberrantExpression_pipeline_Counting_Datasets_R
	2	AberrantExpression_pipeline_Counting_Summary_R
	10	AberrantExpression_pipeline_Counting_countReads_R
	2	AberrantExpression_pipeline_Counting_filterCounts_R
	2	AberrantExpression_pipeline_Counting_mergeCounts_R
	1	AberrantExpression_pipeline_Counting_preprocessGeneAnnotation_R
	1	AberrantExpression_pipeline_OUTRIDER_Datasets_R
	2	AberrantExpression_pipeline_OUTRIDER_Summary_R
	2	AberrantExpression_pipeline_OUTRIDER_results_R
	2	AberrantExpression_pipeline_OUTRIDER_runOutrider_R
	1	AberrantSplicing_Overview_R
	1	AberrantSplicing_pipeline_Counting_00_define_datasets_from_anno_R
	1	AberrantSplicing_pipeline_Counting_01_0_countRNA_init_R
	10	AberrantSplicing_pipeline_Counting_01_1_countRNA_splitReads_samplewise_R
	1	AberrantSplicing_pipeline_Counting_01_2_countRNA_splitReads_merge_R
	10	AberrantSplicing_pipeline_Counting_01_3_countRNA_nonSplitReads_samplewise_R
	1	AberrantSplicing_pipeline_Counting_01_4_countRNA_nonSplitReads_merge_R
	1	AberrantSplicing_pipeline_Counting_01_5_countRNA_collect_R
	1	AberrantSplicing_pipeline_Counting_02_psi_value_calculation_FraseR_R
	1	AberrantSplicing_pipeline_Counting_03_filter_expression_FraseR_R
	1	AberrantSplicing_pipeline_Counting_DatasetsF_R
	1	AberrantSplicing_pipeline_Counting_Summary_R
	1	AberrantSplicing_pipeline_FRASER_04_fit_hyperparameters_FraseR_R
	1	AberrantSplicing_pipeline_FRASER_05_fit_autoencoder_FraseR_R
	1	AberrantSplicing_pipeline_FRASER_06_calculation_stats_AE_FraseR_R
	1	AberrantSplicing_pipeline_FRASER_07_extract_results_FraseR_R
	1	AberrantSplicing_pipeline_FRASER_Datasets_R
	1	AberrantSplicing_pipeline_FRASER_Summary_R
	1	Index
	1	MonoallelicExpression_Overview_R
	1	MonoallelicExpression_pipeline_MAE_Datasets_R
	1	MonoallelicExpression_pipeline_MAE_Results_R
	2	MonoallelicExpression_pipeline_MAE_deseq_mae_R
	1	MonoallelicExpression_pipeline_QC_DNA_RNA_matrix_plot_R
	1	MonoallelicExpression_pipeline_QC_Datasets_R
	1	MonoallelicExpression_pipeline_QC_create_matrix_dna_rna_cor_R
	2	MonoallelicExpression_pipeline_QC_deseq_qc_R
	1	aberrantExpression
	4	aberrantExpression_bamStats
	2	aberrantExpression_mergeBamStats
	1	aberrantSplicing
	1	aberrantSplicing_dependency
	1	all
	1	dependencyGraph
	1	mae
	4	mae_allelicCounts
	4	mae_createSNVs
	1	mae_dependency
	94
Select jobs to execute...

[Thu May 13 17:33:29 2021]
rule mae_dependency:
    output: /gpfs/gsfs10/users/apptest1/Output/html/mae-pipeline_dep.svg
    jobid: 100
...
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
For submitting large jobs through cluster, we will create a cluster.yaml file first. For example, this yaml file (thanks for contribution from Sarah Silverstein) is based on 80 samples about 130M reads per sample (if you have smaller or less samples, please run snakemake at interactive session or scale down accordingly):



```

__default__:
        partition: norm
        threads: 2
        memory: 2G
        time: "2:00:00"

AberrantExpression_pipeline_OUTRIDER_runOutrider_R:
        threads: 4
        memory: 10G
        time: "20:00:00"

AberrantExpression_pipeline_Counting_countReads_R:
        memory: 20G
        time: "6:00:00"

AberrantExpression_pipeline_OUTRIDER_Summary_R:
        threads: 50
        memory: 20G

AberrantExpression_Overview_R:
        memory: 5G


AberrantSplicing_pipeline_Counting_01_1_countRNA_splitReads_samplewise_R:
        memory: 10G
        threads: 4

AberrantSplicing_pipeline_Counting_01_2_countRNA_splitReads_merge_R:
        memory: 10G

AberrantSplicing_pipeline_Counting_01_4_countRNA_nonSplitReads_merge_R:
        memory: 10G

AberrantSplicing_pipeline_Counting_02_psi_value_calculation_FraseR_R:
        memory: 35G
        threads: 20

AberrantSplicing_pipeline_Counting_03_filter_expression_FraseR_R:
        memory: 15G

AberrantSplicing_pipeline_FRASER_04_fit_hyperparameters_FraseR_R:
        memory: 25G
        threads: 10

AberrantSplicing_pipeline_FRASER_05_fit_autoencoder_FraseR_R:
        memory: 50G
        threads: 22
        time: "4:00:00"

AberrantSplicing_pipeline_FRASER_06_calculation_stats_AE_FraseR_R:
        memory: 110G
        threads: 24
        time: "4:00:00"

AberrantSplicing_pipeline_FRASER_07_extract_results_FraseR_R:
        memory: 25G

AberrantSplicing_pipeline_FRASER_Summary_R:
        memory: 15G
        threads: 4

AberrantSplicing_Overview_R:
        memory: 5G


MonoallelicExpression_pipeline_MAE_Results_R:
        memory: 5G

MonoallelicExpression_pipeline_QC_create_matrix_dna_rna_cor_R:
        memory: 200G
        threads: 2
	time: "14:00:00"

mae_allelicCounts:
        threads: 4


```

Create a batch input file (e.g. drop.sh). For example:


hljs.highlightAll();


hljs.highlightAll();

```


#!/bin/bash
set -e
module load drop
cd /data/$USER/drop_test/
snakemake -pr --jobs 10 \
--cluster "sbatch --cpus-per-task={cluster.threads} --mem={cluster.memory} --time={cluster.time}" \
--cluster-config cluster.yaml --latency-wait 120 --max-jobs-per-second 1 \
--max-status-checks-per-second 0.01 all


```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=2 --mem=2g --time=4-00:00:00 drop.sh
```


