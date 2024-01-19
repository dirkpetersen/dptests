

document.querySelector('title').textContent = 'DeepSea on Biowulf';
DeepSea on Biowulf


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



DeepSEA is a deep learning-based algorithmic framework for predicting the chromatin effects of sequence alterations with single nucleotide sensitivity. DeepSEA can accurately predict the epigenetic state of a sequence, including transcription factors binding, DNase I sensitivities and histone marks in multiple cell types, and further utilize this capability to predict the chromatin effects of sequence variants and prioritize regulatory variants.



### References:


* [Zhou, Jian, and Olga G. Troyanskaya. "Predicting effects of noncoding variants with deep learning–based sequence model." *Nature methods* 12.10 (2015): 931-934.](https://www.nature.com/articles/nmeth.3547?source=post_page---------------------------)


Documentation
* [DeepSea Main Site](http://deepsea.princeton.edu/job/analysis/create/)


Important Notes
* Module Name: deepsea (see [the modules page](/apps/modules.html) for more information)
 * Omit the python prefix when using the scripts that make up DeepSea. It is installed within a Singularity container, and the python prefix will be automatically added to ensure that DeepSea uses the appropriate python installed in the container itself. 
 * Example data can be found in /fdb/deepsea/0.94c/example\*. The Gerp, phastCons, and phylop subdirectories that appear in the same location contain reference data and are automatically made available to DeepSea within the container at runtime. There is no need for users to touch these files.
 * The script rundeepsea-insilicomut.py attempts to use all of the CPUs in a node. Therefore users need to allocate nodes exclusively when running this script.
 * The container included here is only useful for classification and does not work for training. Therefore it is a CPU only application. 
 * If you receive an error complaining "No Valid variants found", check that the varients that you've provided are *non-coding variants*. DeepSea is designed to compute functional scores on non-coding variants only. 
 * DeepSea creates temporariy files in the $TMPDIR directory. By default, this directory is set to /tmp. Usersare encouraged to reset this value to /lscratch/$SLURM\_JOB\_ID to avoid filling up /tmp and negatively impacting jobs on a particular node.



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf ~]$ **sinteractive --constraint=x2650 --cpus-per-task=32 --mem=10g --ntasks=1 --exclusive --gres=lscratch:10**
salloc.exe: Pending job allocation 47562568
salloc.exe: job 47562568 queued and waiting for resources
salloc.exe: job 47562568 has been allocated resources
salloc.exe: Granted job allocation 47562568
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn0236 are ready for job
srun: error: x11: no local DISPLAY defined, skipping

[user@cn0236 ~]$ **export TMPDIR=/lscratch/${SLURM\_JOB\_ID}**

[user@cn0236 ~]$ **mkdir /lscratch/${SLURM\_JOB\_ID}/test**

[user@cn0236 ~]$ **cd /lscratch/${SLURM\_JOB\_ID}/test**

[user@cn0236 test]$ **cp /fdb/deepsea/0.94c/example\* .**

[user@cn0236 test]$ **module load deepsea**
[+] Loading deepsea  0.94c  on cn0236
[+] Loading singularity  3.5.2  on cn0236

[user@cn0236 test]$ **rundeepsea.py example.vcf out1**
Successfully copied input to working directory /tmp/tmpIPQUab
Loading required package: BSgenome.Hsapiens.UCSC.hg19
Loading required package: BSgenome
Loading required package: methods
Loading required package: BiocGenerics
Loading required package: parallel

Attaching package: 'BiocGenerics'

The following objects are masked from 'package:parallel':

    clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    clusterExport, clusterMap, parApply, parCapply, parLapply,
    parLapplyLB, parRapply, parSapply, parSapplyLB

The following objects are masked from 'package:stats':

    IQR, mad, xtabs

The following objects are masked from 'package:base':

    Filter, Find, Map, Position, Reduce, anyDuplicated, append,
    as.data.frame, as.vector, cbind, colnames, do.call, duplicated,
    eval, evalq, get, grep, grepl, intersect, is.unsorted, lapply,
    lengths, mapply, match, mget, order, paste, pmax, pmax.int, pmin,
    pmin.int, rank, rbind, rownames, sapply, setdiff, sort, table,
    tapply, union, unique, unlist, unsplit

Loading required package: S4Vectors
Loading required package: stats4
Loading required package: IRanges
Loading required package: GenomeInfoDb
Loading required package: GenomicRanges
Loading required package: Biostrings
Loading required package: XVector
Loading required package: rtracklayer
RangedData with 4 rows and 3 value columns across 3 spaces
     space                 ranges |         ori         mut
  <factor>              <IRanges> | <character> <character>
1     chr1 [109817041, 109818140] |           G           T
2    chr10 [ 23507814,  23508913] |           A           G
3    chr16 [   209160,    210259] |           T           C
4    chr16 [ 52598639,  52599738] |           C           T
                            name
                     <character>
1  [known_CEBP_binding_increase]
2 [known_FOXA2_binding_decrease]
3 [known_GATA1_binding_increase]
4 [known_FOXA1_binding_increase]
Number of valid variants:
4
Number of input variants:
4
Successfully converted to input format
Processing options
Switch to float
1
Processing options
Switch to float
1
Finished running DeepSEA. Now prepare output files...
[W::hts_idx_load2] The index file is older than the data file: ./resources/phastCons/primates_nohuman.tsv.gz.tbi
(4, 4)
Finished creating output file. Now clean up...
Everything done.

[user@cn0236 test]$ **ls out1/**
infile.vcf.out.alt     infile.vcf.out.logfoldchange  infile.vcf.wt1100.fasta.ref.vcf.evoall
infile.vcf.out.diff    infile.vcf.out.ref            infile.vcf.wt1100.fasta.ref.vcf.evo.evalues
infile.vcf.out.evalue  infile.vcf.out.snpclass
infile.vcf.out.funsig  infile.vcf.out.summary

[user@cn0236 test]$ **rundeepsea-insilicomut.py example.vcf out2 469 #this will try to use all of the CPUs on a node**
Successfully copied input to working directory /tmp/tmpCLZ7CN
Loading required package: BSgenome.Hsapiens.UCSC.hg19
Loading required package: BSgenome
Loading required package: methods
Loading required package: BiocGenerics
Loading required package: parallel

Attaching package: 'BiocGenerics'

The following objects are masked from 'package:parallel':

    clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    clusterExport, clusterMap, parApply, parCapply, parLapply,
    parLapplyLB, parRapply, parSapply, parSapplyLB

The following objects are masked from 'package:stats':

    IQR, mad, xtabs

The following objects are masked from 'package:base':

    Filter, Find, Map, Position, Reduce, anyDuplicated, append,
    as.data.frame, as.vector, cbind, colnames, do.call, duplicated,
    eval, evalq, get, grep, grepl, intersect, is.unsorted, lapply,
    lengths, mapply, match, mget, order, paste, pmax, pmax.int, pmin,
    pmin.int, rank, rbind, rownames, sapply, setdiff, sort, table,
    tapply, union, unique, unlist, unsplit

Loading required package: S4Vectors
Loading required package: stats4
Loading required package: IRanges
Loading required package: GenomeInfoDb
Loading required package: GenomicRanges
Loading required package: Biostrings
Loading required package: XVector
Loading required package: rtracklayer
RangedData with 4 rows and 3 value columns across 3 spaces
     space                 ranges |         ori         mut
  <factor>              <IRanges> | <character> <character>
1     chr1 [109817041, 109818140] |           G           T
2    chr10 [ 23507814,  23508913] |           A           G
3    chr16 [   209160,    210259] |           T           C
4    chr16 [ 52598639,  52599738] |           C           T
                            name
                     <character>
1  [known_CEBP_binding_increase]
2 [known_FOXA2_binding_decrease]
3 [known_GATA1_binding_increase]
4 [known_FOXA1_binding_increase]
Number of valid variants:
4
Number of input variants:
4
Successfully converted to input format
Processing options
Switch to float
1
Processing options
Switch to float
1
1025
2049
3073
4097
5121
6145
7169
8193
9217
10241
11265
12289
13313
14337
15361
16385
17409
18433
19457
20481
21505
22529
23553
Finished running DeepSEA. Now prepare output files...
/opt/conda/lib/python2.7/site-packages/matplotlib/collections.py:590: FutureWarning: elementwise comparison failed; returning scalar instead, but in the future will perform elementwise comparison
  if self._edgecolors == str('face'):
Finished creating output file. Now clean up...
Everything done.

[user@cn0236 test]$ **ls out2/**
colorbar.png  log2foldchange_profile.csv  preview.png  vis.png

[user@cn0236 test]$ **cp -r out\* /data/${USER}**

[user@cn0236 test]$ **exit**
exit
salloc.exe: Relinquishing job allocation 47562568

[user@biowulf ~]$ 

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. deepsea.sh). For example:



```

#!/bin/bash
set -e
module load deepsea
export TMPDIR=/lscratch/$SLURM_JOB_ID
rundeepsea-insilicomut.py example.vcf out

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command. As in the example above, be sure to allocate a node exclusively. For instance, you could submit this job using the following command



```
sbatch --constraint=x2650 --cpus-per-task=32 --mem=10g --ntasks=1 --exclusive --gres=lscratch:10 deepsea.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. deepsea.swarm). For example:



```

export TMPDIR=/lscratch/$SLURM_JOB_ID; rundeepsea-insilicomut.py example1.vcf out1
export TMPDIR=/lscratch/$SLURM_JOB_ID; rundeepsea-insilicomut.py example2.vcf out2
export TMPDIR=/lscratch/$SLURM_JOB_ID; rundeepsea-insilicomut.py example3.vcf out3
export TMPDIR=/lscratch/$SLURM_JOB_ID; rundeepsea-insilicomut.py example4.vcf out4

```

Submit this job using the [swarm](/apps/swarm.html) command. If you are using the rundeepsea-insilicomut.py script, be sure to allocate the nodes exclusively. For example:



```
swarm -f deepsea.swarm -g 10 -t 32 --module deepsea --sbatch "--constraint=x2650 --ntasks=1 --exclusive --gres=lscratch:10"
```








