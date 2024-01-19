

document.querySelector('title').textContent = 'biobakery\_workflows: workflows and tasks for executing common microbial community analyses. ';
**biobakery\_workflows: workflows and tasks for executing common microbial community analyses.** 


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



bioBakery is a meta’omic analysis environment and collection of individual software
tools with the capacity to process raw shotgun sequencing data into actionable microbial community
feature profiles, summary reports, and publication-ready figures. 
It includes a collection of preconfigured analysis modules also joined into workflows for reproducibility.
Each individual module has been developed to perform a particular task, e.g. quantitative taxonomic
profiling or statistical analysis.



References:
-----------


* Lauren J McIver, Galeb Abu-Ali, Eric A Franzosa, Randall Schwager, Xochitl C Morgan, 
 Levi Waldron, Nicola Segata, Curtis Huttenhower   

 *bioBakery: a meta’omic analysis environment*   

[Bioinformatics](https://academic.oup.com/bioinformatics/article/34/7/1235/4673198?login=true)  
 Volume 34, Issue 7, 01 April 2018, Pages 1235–1237


Documentation
* [biobakery\_workflows Github page](https://github.com/biobakery/biobakery_workflows)
* [bioBakery workflows home page](https://huttenhower.sph.harvard.edu/biobakery_workflows)


Important Notes
* Module Name: biobakery\_workflows (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set
	+ **BW\_HOME**  installation directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
  


```

[user@biowulf]$ **sinteractive** 
[user@cn0861 ~]$ **module load biobakery\_workflows** 
[+] Loading samtools 1.15.1  ...
[+] Loading bowtie  2-2.4.5
[+] Loading trimmomatic  0.39  on cn4276
[+] Loading singularity  3.10.5  on cn4276
[+] Loading metaphlan  4.0.3
[+] Loading biobakery_workflows 3.1  ...

```

All workflows follow the general command format:

```

biobakery_workflows $WORKFLOW --input $INPUT --output $OUTPUT

```

where $WORKFLOW is one of: 

```

16s, 16s_vis, isolate_assembly, wmgx, wmgx_vis, wmgx_wmtx, wmgx_wmtx_vis

```

The basic usage of the biobakery\_workflows executable is as follows:

```

[user@cn0861 ~]$ **biobakery\_workflows -h**
usage: biobakery_workflows [-h] [--version]
                           {16s,16s_vis,isolate_assembly,wmgx,wmgx_vis,wmgx_wmtx,wmgx_wmtx_vis}

bioBakery workflows: A collection of AnADAMA2 workflows

positional arguments:
  {16s,16s_vis,isolate_assembly,wmgx,wmgx_vis,wmgx_wmtx,wmgx_wmtx_vis}
                        workflow to run

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit

```

To see the usage of biobakery\_workflows for a particular workflow, say, "wmgx", enter the command:

```

[user@cn0861 ~]$ **biobakery\_workflows wmgx -h**
usage: wmgx.py [-h] [--version]
               [--input-extension {fastq.gz,fastq,fq.gz,fq,fasta,fasta.gz,fastq.bz2,fq.bz2}]
               [--barcode-file BARCODE_FILE]
               [--dual-barcode-file DUAL_BARCODE_FILE]
               [--index-identifier INDEX_IDENTIFIER]
               [--min-pred-qc-score MIN_PRED_QC_SCORE] [--threads THREADS]
               [--pair-identifier PAIR_IDENTIFIER] [--interleaved]
               [--bypass-quality-control]
               [--contaminate-databases CONTAMINATE_DATABASES]
               [--qc-options QC_OPTIONS]
               [--functional-profiling-options FUNCTIONAL_PROFILING_OPTIONS]
               [--remove-intermediate-output] [--bypass-functional-profiling]
               [--bypass-strain-profiling] [--run-strain-gene-profiling]
               [--bypass-taxonomic-profiling] [--run-assembly]
               [--strain-profiling-options STRAIN_PROFILING_OPTIONS]
               [--taxonomic-profiling-options TAXONOMIC_PROFILING_OPTIONS]
               [--max-strains MAX_STRAINS] [--strain-list STRAIN_LIST]
               [--assembly-options ASSEMBLY_OPTIONS] -o OUTPUT [-i INPUT]
               [--config CONFIG] [--local-jobs JOBS] [--grid-jobs GRID_JOBS]
               [--grid GRID] [--grid-partition GRID_PARTITION]
               [--grid-benchmark {on,off}] [--grid-options GRID_OPTIONS]
               [--grid-environment GRID_ENVIRONMENT]
               [--grid-scratch GRID_SCRATCH] [--dry-run] [--skip-nothing]
               [--quit-early] [--until-task UNTIL_TASK]
               [--exclude-task EXCLUDE_TASK] [--target TARGET]
               [--exclude-target EXCLUDE_TARGET]
               [--log-level {DEBUG,INFO,WARNING,ERROR,CRITICAL}]

A workflow for whole metagenome shotgun sequences

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
  --input-extension {fastq.gz,fastq,fq.gz,fq,fasta,fasta.gz,fastq.bz2,fq.bz2}
                        the input file extension
                        [default: fastq.gz]
  --barcode-file BARCODE_FILE
                        the barcode file
                        [default: ]
  --dual-barcode-file DUAL_BARCODE_FILE
                        the string to identify the dual barcode file
                        [default: ]
  --index-identifier INDEX_IDENTIFIER
                        the string to identify the index files
                        [default: _I1_001]
  --min-pred-qc-score MIN_PRED_QC_SCORE
                        the min phred quality score to use for demultiplexing
                        [default: 2]
  --threads THREADS     number of threads/cores for each task to use
                        [default: 1]
  --pair-identifier PAIR_IDENTIFIER
                        the string to identify the first file in a pair
                        [default: .R1]
  --interleaved         indicates whether or not sequence files are interleaved
                        [default: False]
  --bypass-quality-control
                        do not run the quality control tasks
  --contaminate-databases CONTAMINATE_DATABASES
                        the path (or comma-delimited paths) to the contaminate
                        reference databases for QC
                        [default: /fdb/biobakery_workflows_databases/kneaddata_db_human_genome]
  --qc-options QC_OPTIONS
                        additional options when running the QC step
                        [default: ]
  --functional-profiling-options FUNCTIONAL_PROFILING_OPTIONS
                        additional options when running the functional profiling step
                        [default: ]
  --remove-intermediate-output
                        remove intermediate output files
  --bypass-functional-profiling
                        do not run the functional profiling tasks
  --bypass-strain-profiling
                        do not run the strain profiling tasks (StrainPhlAn)
  --run-strain-gene-profiling
                        run the gene-based strain profiling tasks (PanPhlAn)
  --bypass-taxonomic-profiling
                        do not run the taxonomic profiling tasks (a tsv profile for each sequence file must be included in the input folder using the same sample name)
  --run-assembly        run the assembly and annotation tasks
  --strain-profiling-options STRAIN_PROFILING_OPTIONS
                        additional options when running the strain profiling step
                        [default: ]
  --taxonomic-profiling-options TAXONOMIC_PROFILING_OPTIONS
                        additional options when running the taxonomic profiling step
                        [default: ]
  --max-strains MAX_STRAINS
                        the max number of strains to profile
                        [default: 20]
  --strain-list STRAIN_LIST
                        input file with list of strains to profile
                        [default: ]
  --assembly-options ASSEMBLY_OPTIONS
                        additional options when running the assembly step
                        [default: ]
  -o OUTPUT, --output OUTPUT
                        Write output to this directory
  -i INPUT, --input INPUT
                        Find inputs in this directory
                        [default: /gpfs/gsfs7/users/$USER/biobakery_workflows]
  --config CONFIG       Find workflow configuration in this folder
                        [default: only use command line options]
  --local-jobs JOBS     Number of tasks to execute in parallel locally
                        [default: 1]
  --grid-jobs GRID_JOBS
                        Number of tasks to execute in parallel on the grid
                        [default: 0]
  --grid GRID           Run gridable tasks on this grid type
                        [default: slurm]
  --grid-partition GRID_PARTITION
                        Partition/queue used for gridable tasks.
                        Provide a single partition or a comma-delimited list
                        of short/long partitions with a cutoff.
                        [default: serial_requeue,shared,240]
  --grid-benchmark {on,off}
                        Benchmark gridable tasks
                        [default: on]
  --grid-options GRID_OPTIONS
                        Grid specific options that will be applied to each grid task
  --grid-environment GRID_ENVIRONMENT
                        Commands that will be run before each grid task to set up environment
  --grid-scratch GRID_SCRATCH
                        The folder to write intermediate scratch files for grid jobs
  --dry-run             Print tasks to be run but don't execute their actions
  --skip-nothing        Run all tasks. Rerun tasks that have already been run.
  --quit-early          Stop if a task fails. By default,
                        all tasks (except sub-tasks of failed tasks) will run.
  --until-task UNTIL_TASK
                        Stop after running this task. Use task name or number.
  --exclude-task EXCLUDE_TASK
                        Don't run these tasks. Add multiple times to append.
  --target TARGET       Only run tasks that generate these targets.
                        Add multiple times to append.
                        Patterns with ? and * are allowed.
  --exclude-target EXCLUDE_TARGET
                        Don't run tasks that generate these targets.
                        Add multiple times to append.
                        Patterns with ? and * are allowed.
  --log-level {DEBUG,INFO,WARNING,ERROR,CRITICAL}
                        Set the level of output for the log
                        [default: INFO]

```

This output indicates, in particular, that the biobakery executable can be run with multiple threads (option --threads THREADS) if an appropriate number of cpus is allocated for the interactive session, which should speed up your processing.   
   

Besides, there's as number of other relevant executables available:

```

[user@cn0861 ~]$ **ls $BW\_BIN**
16s.py                                        merge_and_rename_fastq.py
16s_vis.py                                    merge_fastq.py
anadama2_add_files_to_database.py             merge_paired_ends.R
annotate_genome.py                            phylogeny.R
assign_taxonomy.R                             pull_out_reads_by_species_metaphlan2_results.py
biobakery_workflows                           python
burst_workflow.py                             R
check_fastq_format.py                         remove_if_exists.py
const_seq_table.R                             rename_data_products.py
count_features.py                             rename_fastq_files.py
create_fasta_per_taxonomy_from_alignments.py  rename_files_to_sample_ids.py
create_otu_tables_from_alignments.py          reverse_compliment_barcodes.py
create_subsampled_demos.py                    rna_dna_norm.py
dada2_version.R                               shell
demultiplex_split_index.py                    strainphlan_ggtree_vis.R
extract_orphan_reads.py                       strainphlan_ordination_vis.R
filter_and_trim.R                             trim_taxonomy.py
generate_dual_barcode.py                      wmgx.py
get_counts_from_humann2_logs.py               wmgx_vis.py
identify_primers.R                            wmgx_wmtx.py
isolate_assembly.py.py                        wmgx_wmtx_vis.py
learn_error_rates.R

```

In order to run biobakery\_workflows on sample fastq or fastq.gz data (single-end or paired-end), enter the following commands: 

```

[user@cn0861 ~]$ **biobakery\_workflows wmgx --input $BW\_EXAMPLES/wmgx/single/ --output /data/$USER/workflow\_output**
(Apr 28 09:03:50) [ 0/52 -   0.00%] **Ready    ** Task  3: kneaddata____demo2
(Apr 28 09:03:50) [ 0/52 -   0.00%] **Started  ** Task  3: kneaddata____demo2
(Apr 28 09:03:58) [ 1/52 -   1.92%] **Completed** Task  3: kneaddata____demo2
(Apr 28 09:03:58) [ 1/52 -   1.92%] **Ready    ** Task  8: metaphlan____demo2
(Apr 28 09:03:58) [ 1/52 -   1.92%] **Started  ** Task  8: metaphlan____demo2
(Apr 28 09:05:41) [ 2/52 -   3.85%] **Completed** Task  8: metaphlan____demo2
(Apr 28 09:05:41) [ 2/52 -   3.85%] **Ready    ** Task 34: strainphlan_sample2markers____demo2
(Apr 28 09:05:41) [ 2/52 -   3.85%] **Started  ** Task 34: strainphlan_sample2markers____demo2
(Apr 28 09:06:18) [ 3/52 -   5.77%] **Completed** Task 34: strainphlan_sample2markers____demo2
(Apr 28 09:06:18) [ 3/52 -   5.77%] **Ready    ** Task 13: humann____demo2
(Apr 28 09:06:18) [ 3/52 -   5.77%] **Started  ** Task 13: humann____demo2
...
(Apr 28 09:12:25) [49/52 -  94.23%] **Completed** Task 24: humann_renorm_pathways_relab____demo1
(Apr 28 09:12:25) [49/52 -  94.23%] **Ready    ** Task 28: humann_join_tables_pathways_relab
(Apr 28 09:12:25) [49/52 -  94.23%] **Started  ** Task 28: humann_join_tables_pathways_relab
(Apr 28 09:12:25) [50/52 -  96.15%] **Completed** Task 28: humann_join_tables_pathways_relab
(Apr 28 09:12:25) [50/52 -  96.15%] **Ready    ** Task 31: humann_count_features_pathways
(Apr 28 09:12:25) [50/52 -  96.15%] **Started  ** Task 31: humann_count_features_pathways
(Apr 28 09:12:25) [51/52 -  98.08%] **Completed** Task 31: humann_count_features_pathways
(Apr 28 09:12:25) [51/52 -  98.08%] **Ready    ** Task 32: humann_merge_feature_counts
(Apr 28 09:12:25) [51/52 -  98.08%] **Started  ** Task 32: humann_merge_feature_counts
(Apr 28 09:12:25) [52/52 - 100.00%] **Completed** Task 32: humann_merge_feature_counts
Run Finished
[user@cn0861 ~]$ **tree /data/user/workflow\_output**
/data/user/workflow_output
├── anadama.log
├── humann
│   ├── counts
│   │   ├── humann_ecs_relab_counts.tsv
│   │   ├── humann_feature_counts.tsv
│   │   ├── humann_genefamilies_relab_counts.tsv
│   │   ├── humann_pathabundance_relab_counts.tsv
│   │   └── humann_read_and_species_count_table.tsv
│   ├── main
│   │   ├── demo1_genefamilies.tsv
│   │   ├── demo1_humann_temp
│   │   │   ├── demo1_bowtie2_aligned.sam
│   │   │   ├── demo1_bowtie2_aligned.tsv
│   │   │   ├── demo1_bowtie2_index.1.bt2
│   │   │   ├── demo1_bowtie2_index.2.bt2
│   │   │   ├── demo1_bowtie2_index.3.bt2
│   │   │   ├── demo1_bowtie2_index.4.bt2
│   │   │   ├── demo1_bowtie2_index.rev.1.bt2
│   │   │   ├── demo1_bowtie2_index.rev.2.bt2
│   │   │   ├── demo1_bowtie2_unaligned.fa
│   │   │   ├── demo1_custom_chocophlan_database.ffn
│   │   │   ├── demo1_diamond_aligned.tsv
│   │   │   └── demo1_diamond_unaligned.fa
│   │   ├── demo1.log
│   │   ├── demo1_pathabundance.tsv
│   │   ├── demo1_pathcoverage.tsv
│   │   ├── demo2_genefamilies.tsv
│   │   ├── demo2_humann_temp
│   │   │   ├── demo2_bowtie2_aligned.sam
│   │   │   ├── demo2_bowtie2_aligned.tsv
│   │   │   ├── demo2_bowtie2_index.1.bt2
│   │   │   ├── demo2_bowtie2_index.2.bt2
│   │   │   ├── demo2_bowtie2_index.3.bt2
│   │   │   ├── demo2_bowtie2_index.4.bt2
│   │   │   ├── demo2_bowtie2_index.rev.1.bt2
│   │   │   ├── demo2_bowtie2_index.rev.2.bt2
│   │   │   ├── demo2_bowtie2_unaligned.fa
│   │   │   ├── demo2_custom_chocophlan_database.ffn
│   │   │   ├── demo2_diamond_aligned.tsv
│   │   │   └── demo2_diamond_unaligned.fa
│   │   ├── demo2.log
│   │   ├── demo2_pathabundance.tsv
│   │   └── demo2_pathcoverage.tsv
│   ├── merged
│   │   ├── ecs_relab.tsv
│   │   ├── ecs.tsv
│   │   ├── genefamilies_relab.tsv
│   │   ├── genefamilies.tsv
│   │   ├── pathabundance_relab.tsv
│   │   └── pathabundance.tsv
│   ├── regrouped
│   │   ├── demo1_ecs.tsv
│   │   └── demo2_ecs.tsv
│   └── relab
│       ├── ecs
│       │   ├── demo1_ecs_relab.tsv
│       │   └── demo2_ecs_relab.tsv
│       ├── genes
│       │   ├── demo1_genefamilies_relab.tsv
│       │   └── demo2_genefamilies_relab.tsv
│       └── pathways
│           ├── demo1_pathabundance_relab.tsv
│           └── demo2_pathabundance_relab.tsv
├── kneaddata
│   ├── main
│   │   ├── demo1.fastq
│   │   ├── demo1_Homo_sapiens_bowtie2_contam.fastq
│   │   ├── demo1.log
│   │   ├── demo1.trimmed.fastq
│   │   ├── demo2.fastq
│   │   ├── demo2_Homo_sapiens_bowtie2_contam.fastq
│   │   ├── demo2.log
│   │   └── demo2.trimmed.fastq
│   └── merged
│       └── kneaddata_read_count_table.tsv
├── metaphlan
│   ├── main
│   │   ├── demo1_bowtie2.sam
│   │   ├── demo1_taxonomic_profile.tsv
│   │   ├── demo2_bowtie2.sam
│   │   └── demo2_taxonomic_profile.tsv
│   └── merged
│       ├── metaphlan_species_counts_table.tsv
│       └── metaphlan_taxonomic_profiles.tsv
└── strainphlan
    ├── 0_clade.log
    ├── 0_clade.tree
    ├── 10_clade.log
    ├── 10_clade.tree
    ├── 11_clade.log
   ...
    ├── 9_clade.log
    ├── 9_clade.tree
    ├── clades_list_order_by_average_abundance.txt
    ├── clades_list.txt
    ├── demo1_bowtie2
    │   └── demo1_bowtie2.pkl
    └── demo2_bowtie2
        └── demo2_bowtie2.pkl
[user@cn0861 ~]$ **biobakery\_workflows wmgx --input $BW\_EXAMPLES/wmgx/paired/ --output /data/$USER/workflow\_output**
(Apr 28 09:14:52) [ 0/52 -   0.00%] **Ready    ** Task  4: kneaddata____demo2
(Apr 28 09:14:52) [ 0/52 -   0.00%] **Started  ** Task  4: kneaddata____demo2
(Apr 28 09:15:05) [ 1/52 -   1.92%] **Completed** Task  4: kneaddata____demo2
(Apr 28 09:15:05) [ 1/52 -   1.92%] **Ready    ** Task 10: metaphlan____demo2
(Apr 28 09:15:05) [ 1/52 -   1.92%] **Started  ** Task 10: metaphlan____demo2
(Apr 28 09:16:49) [ 2/52 -   3.85%] **Completed** Task 10: metaphlan____demo2
(Apr 28 09:16:49) [ 2/52 -   3.85%] **Ready    ** Task 36: strainphlan_sample2markers____demo2
(Apr 28 09:16:49) [ 2/52 -   3.85%] **Started  ** Task 36: strainphlan_sample2markers____demo2
...
(Apr 28 09:27:03) [46/52 -  88.46%] **Completed** Task 22: humann_renorm_genes_relab____demo1
(Apr 28 09:27:03) [46/52 -  88.46%] **Ready    ** Task 28: humann_join_tables_genes_relab
(Apr 28 09:27:03) [46/52 -  88.46%] **Started  ** Task 28: humann_join_tables_genes_relab
(Apr 28 09:27:03) [47/52 -  90.38%] **Completed** Task 28: humann_join_tables_genes_relab
(Apr 28 09:27:03) [47/52 -  90.38%] **Ready    ** Task 31: humann_count_features_genes
(Apr 28 09:27:03) [47/52 -  90.38%] **Started  ** Task 31: humann_count_features_genes
(Apr 28 09:27:03) [48/52 -  92.31%] **Completed** Task 31: humann_count_features_genes
(Apr 28 09:27:03) [48/52 -  92.31%] **Ready    ** Task 26: humann_renorm_pathways_relab____demo1
(Apr 28 09:27:03) [48/52 -  92.31%] **Started  ** Task 26: humann_renorm_pathways_relab____demo1
(Apr 28 09:27:03) [49/52 -  94.23%] **Completed** Task 26: humann_renorm_pathways_relab____demo1
(Apr 28 09:27:03) [49/52 -  94.23%] **Ready    ** Task 30: humann_join_tables_pathways_relab
(Apr 28 09:27:03) [49/52 -  94.23%] **Started  ** Task 30: humann_join_tables_pathways_relab
(Apr 28 09:27:03) [50/52 -  96.15%] **Completed** Task 30: humann_join_tables_pathways_relab
(Apr 28 09:27:03) [50/52 -  96.15%] **Ready    ** Task 33: humann_count_features_pathways
(Apr 28 09:27:03) [50/52 -  96.15%] **Started  ** Task 33: humann_count_features_pathways
(Apr 28 09:27:03) [51/52 -  98.08%] **Completed** Task 33: humann_count_features_pathways
(Apr 28 09:27:03) [51/52 -  98.08%] **Ready    ** Task 34: humann_merge_feature_counts
(Apr 28 09:27:03) [51/52 -  98.08%] **Started  ** Task 34: humann_merge_feature_counts
(Apr 28 09:27:03) [52/52 - 100.00%] **Completed** Task 34: humann_merge_feature_counts
Run Finished



```


