

document.querySelector('title').textContent = 'VirSorter2: a multi-classifier, expert-guided approach to detect diverse DNA and RNA viruses';
**VirSorter2: a multi-classifier, expert-guided approach to detect diverse DNA and RNA viruses**


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



VirSorter2 is a DNA and RNA virus identification tool. It leverages genome-informed
database advances across a collection of customized automatic classifiers 
to improve the accuracy and range of virus sequence detection.



Documentation
* [VirSorter2 Github / tutorial page](https://github.com/jiarong/VirSorter2)


Important Notes
* Module Name: VirSorter2 (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set
	+ **VS2\_HOME**  installation directory
	+ **VS2\_BIN**       executable directory
	+ **VS2\_DATA**    sample data folder



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. 
  
Sample session on a GPU node:



```

[user@biowulf ~]$ **sinteractive --mem=100g --cpus-per-task=4 --gres=lscratch:50**
[user@cn2379 ~]$ **module load virsorter2**
[+] Loading singularity  3.8.2  on cn0853
[+] Loading virsorter2  2.2.3
[user@cn2379 ~]$ **mkdir -p /data/$USER/VirSorter2 && cd /data/$USER/VirSorter2**

```

Configuring VirSorter2 (need to be done only once):

```
 
[user@cn2379 ~]$ **wget https://osf.io/v46sc/download -O db.tgz** 
[user@cn2379 ~]$ **gunzip -c db.tgz | tar -xvf -**
[user@cn2379 ~]$ **vs2 config --init-source --db-dir=$PWD/db**
[2021-09-30 06:23 INFO] VirSorter 2.2.3
[2021-09-30 06:23 INFO] /VirSorter/bin/virsorter config --init-source --db-dir=/data/user/db
[2021-09-30 06:23 INFO] Attention: can not write template-config.yaml in source directory:
/VirSorter/lib/python3.7/site-packages/virsorter
makeing a copy to user home direcotry:
/home/user/.virsorter/template-config.yaml

[2021-09-30 06:23 INFO] Using {template} as config template
[2021-09-30 06:23 INFO] saving /data/user/VirSorter2/db as DBDIR to config file /home/user/.virsorter/template-config.yaml

```

Running VirSorter2 on sample data:

```

[user@cn2379 ~]$ **cp $VS2\_DATA/\* .**
[user@cn2379 ~]$ **vs2 run -w test.out -i test-for-sop.fa --min-length 1500 -j 4 all** 
[2021-09-30 17:03 INFO] VirSorter 2.1
[2021-09-30 17:03 INFO] /miniconda/bin/virsorter run -w test.out -i test-for-sop.fa --min-length 1500 -j 4 all
[2021-10-04 06:24 INFO] VirSorter 2.2.3
[2021-10-04 06:24 INFO] /opt/conda/envs/vs2/bin/virsorter config --init-source --db-dir=./db
[2021-10-04 06:24 INFO] Attention: can not write template-config.yaml in source directory:
/opt/conda/envs/vs2/lib/python3.8/site-packages/virsorter
makeing a copy to user home direcotry:
/home/user/.virsorter/template-config.yaml

[2021-10-04 06:24 INFO] Using {template} as config template
[2021-10-04 06:24 INFO] saving /gs7/users/user/VirSorter2/db as DBDIR to config file /home/user/.virsorter/template-config.yaml
[user@cn0864 VirSorter2]$ ls db/conda_envs
[userga@cn0864 VirSorter2]$  vs2 run -w test.out -i test-for-sop.fa --min-length 1500 -j 4 all
[2021-10-04 06:24 INFO] VirSorter 2.2.3
[2021-10-04 06:24 INFO] /opt/conda/envs/vs2/bin/virsorter run -w test.out -i test-for-sop.fa --min-length 1500 -j 4 all
[2021-10-04 06:24 INFO] Using /home/user/.virsorter/template-config.yaml as config template
[2021-10-04 06:24 INFO] conig file written to /gs7/users/user/VirSorter2/test.out/config.yaml

[2021-10-04 06:24 INFO] Executing: snakemake --snakefile /opt/conda/envs/vs2/lib/python3.8/site-packages/virsorter/Snakefile --directory /gs7/users/user/VirSorter2/test.out --jobs 4 --configfile /gs7/users/user/VirSorter2/test.out/config.yaml --latency-wait 600 --rerun-incomplete --nolock  --conda-frontend mamba --conda-prefix /gs7/users/user/VirSorter2/db/conda_envs --use-conda    --quiet  all
Job counts:
        count   jobs
        1       all
        1       check_point_for_reclassify
        1       circular_linear_split
        1       classify
        2       classify_by_group
        2       classify_full_and_part_by_group
        1       combine_linear_circular
        2       combine_linear_circular_by_group
        1       extract_feature
        1       extract_provirus_seqs
        1       finalize
        1       gff_feature
        2       gff_feature_by_group
        2       hmm_features_by_group
        1       hmm_sort_to_best_hit_taxon
        2       hmm_sort_to_best_hit_taxon_by_group
        1       merge_classification
        1       merge_full_and_part_classification
        2       merge_hmm_gff_features_by_group
        2       merge_provirus_call_by_group_by_split
        1       merge_provirus_call_from_groups
        5       merge_split_hmmtbl
        10      merge_split_hmmtbl_by_group
        10      merge_split_hmmtbl_by_group_tmp
        1       pick_viral_fullseq
        1       preprocess
        1       split_faa
        2       split_faa_by_group
        2       split_gff_by_group
        61
[2021-10-04 06:26 INFO] # of seqs < 1500 bp and removed: 0
[2021-10-04 06:26 INFO] # of circular seqs: 0
[2021-10-04 06:26 INFO] # of linear seqs  : 7
[2021-10-04 06:26 INFO] No circular seqs found in contig file
[2021-10-04 06:26 INFO] Finish spliting linear contig file with common rbs
[2021-10-04 06:26 INFO] Step 1 - preprocess finished.
[2021-10-04 06:32 INFO] Step 2 - extract-feature finished.
[2021-10-04 06:33 INFO]
            ====> VirSorter run (provirus mode) finished.
            # of full    seqs (>=2 genes) as viral:     6
            # of partial seqs (>=2 genes) as viral:     1
            # of short   seqs (< 2 genes) as viral:     0

            Useful output files:
                final-viral-score.tsv       ==> score table
                final-viral-combined.fa     ==> all viral seqs
                final-viral-boundary.tsv    ==> table with boundary info


                Suffix is added to seq names in final-viral-combined.fa:
                full    seqs (>=2 genes) as viral:      ||full
                partial seqs (>=2 genes) as viral:      ||partial
                short   seqs (< 2 genes) as viral:      ||lt2gene


            NOTES:
            Users can further screen the results based on the following
                columns in final-viral-score.tsv:
                - contig length (length)
                - hallmark gene count (hallmark)
                - viral gene % (viral)
                - cellular gene % (cellular)
            The group field in final-viral-score.tsv should NOT be used
                as reliable taxonomy info

            <====

[2021-10-04 06:33 INFO] Step 3 - classify finished.

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. virsorter2.sh). For example:



```

#!/bin/bash
set -e
module load virsorter2
cp $VS2_DATA/* .  
vs2 run -w test.out -i test-for-sop.fa --min-length 1500 -j 4 all

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
**sbatch virsorter2.sh**
```







