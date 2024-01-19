

document.querySelector('title').textContent = 'MINTIE: identifying novel, rare transcriptional variants in cancer RNA-seq data ';
**MINTIE: identifying novel, rare transcriptional variants in cancer RNA-seq data** 


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



MINTIE is a tool for identifying novel, rare transcriptional variants in cancer RNA-seq data. MINTIE detects 
gene fusions, transcribed structural variants, novel splice variants and complex variants, and annotates all novel transcriptional variants.



### References:


* Marek Cmero, Breon Schmidt, Ian J. Majewski, Paul G. Ekert, Alicia Oshlack & Nadia M. Davidson   

*MINTIE: identifying novel structural and splice variants in transcriptomes using RNAseq data*   

[Genome Biology](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02507-8), (2021) 22:296; https://doi.org/10.1186/s13059-021-02507-8


Documentation
* [MINTIE Github page](https://github.com/Oshlack/MINTIE)
* [MINTIE Wiki page](https://github.com/Oshlack/MINTIE/wiki)


Important Notes
* Module Name: MINTIE (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set
	+ **MINTIE\_HOME**  installation directory
	+ **MINTIE\_BIN**       executable directory
	+ **MINTIE\_DATA**  sample data directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=16g -c 16 --gres=lscratch:20**
[user@cn3335 ~]$**module load MINTIE** 
[+] Loading java 12.0.1  ...
[+] Loading MINTIE 0.4.2 ...
[user@cn3335 ~]$ **mkdir -p /data/$USER/MINTIE && cd /data/$USER/MINTIE** 

```

Produce testing data:

```

[user@cn3335 ~]$**mintie -t**

```

and run mintie on the these data: 

```

[user@cn3335 ~]$ **mintie -w -p test\_params.txt cases/\*.fastq.gz controls/\*.fastq.gz**   
????????????????????????????????????????????????????????????????????????????????????????????????????
|                              Starting Pipeline at 2022-12-02 13:05                               |
????????????????????????????????????????????????????????????????????????????????????????????????????

================================ Stage fastq_dedupe (allvars-case) =================================
...
==================================== Stage trim (allvars-case) =====================================
...
================================== Stage assemble (allvars-case) ===================================
...
============================= Stage create_salmon_index (allvars-case) =============================
...
================================= Stage run_salmon (allvars-case) ==================================
...
================================ Stage run_salmon (allvars-control) ================================
...
=========================== Stage create_ec_count_matrix (allvars-case) ============================
...
=================================== Stage run_de (allvars-case) ====================================
...
========================== Stage filter_on_significant_ecs (allvars-case) ==========================
...
======================== Stage align_contigs_against_genome (allvars-case) =========================
...
============================= Stage sort_and_index_bam (allvars-case) ==============================
...
============================== Stage annotate_contigs (allvars-case) ===============================
...
=============================== Stage refine_contigs (allvars-case) ================================
...
================================ Stage calculate_VAF (allvars-case) ================================
...
================================ Stage post_process (allvars-case) =================================
index file allvars-case/novel_contigs.fasta.fai not found, generating...

======================================== Pipeline Succeeded ========================================
20:30:00 MSG:  Finished at Tue Jul 25 20:30:00 EDT 2023
20:30:00 MSG:  Outputs are:
                allvars-case/vaf_estimates.txt
                allvars-case/allvars-case_results.tsv
                allvars-case/novel_contigs.bam
                allvars-case/novel_contigs.vcf
                allvars-case/all_fasta_index/allvars-case.fasta (pre-existing)
                ... 4 more ...
[user@cn3335 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```





