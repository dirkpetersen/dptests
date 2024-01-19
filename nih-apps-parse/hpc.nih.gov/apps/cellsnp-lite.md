

document.querySelector('title').textContent = ' cellSNP-lite on Biowulf';
 cellSNP-lite on Biowulf


|  |
| --- |
| 
Quick Links
[Interactive Jobs](#int)
[Batch job on Biowulf](#batch)
[Using Swarm](#swarm)
[Documentation](#doc)
 |



cellsnp-lite aims to pileup the expressed alleles in single-cell or bulk RNA-seq data, which can be directly used for donor deconvolution in multiplexed single-cell RNA-seq data, particularly with vireo, which assigns cells to donors and detects doublets, even without genotyping reference.

cellsnp-lite heavily depends on htslib. This program should give very similar results as samtools/bcftools mpileup. Also, there are two major differences comparing to bcftools mpileup:

cellsnp-lite can now pileup a list of positions, with directly splitting into a list of cell barcodes, e.g., for 10x genome. With bcftools, you may need to manipulate the RG tag in the bam file if you want to divide reads into cell barcode groups.
cellsnp-lite uses simple filtering for outputting SNPs, i.e., total UMIs or counts and minor alleles fractions. The idea here is to keep most information of SNPs and the downstream statistical model can take the full use of it.
cellsnp-lite is the C version of cellSNP, which is implemented in Python. Compared to cellSNP, cellsnp-lite is basically more efficient with higher speed and less memory usage.

Documentation
<https://github.com/single-cell-genetics/cellsnp-lite>


Important Notes
* Module Name: cellsnp-lite (see [the modules page](/apps/modules.html) for more information)
* Example files are under /usr/local/apps/cellsnp-lite/test
* Mutithreaded, use -p $SLURM\_CPUS\_PER\_TASK to specify number of cpus.


Submitting an interactive job

Allocate an interactive session and run the interactive job there.



```
[biowulf]$ **sinteractive --mem=5g --cpus-per-task=4**
salloc.exe: Granted job allocation 789523
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn0135 are ready for job

[cn0135]$ **cd /data/$USER/**

[cn0135]$ **module load cellsnp-lite**

[cn0135]$ **cp -r /usr/local/apps/cellsnp-lite/test .**

[cn0135]$ **cd test**

[cn0135]$ **bash ./test\_10x.sh**

[cn0135]$ **exit**
salloc.exe: Job allocation 789523 has been revoked.
[biowulf]$
```


Submitting a single batch job
1. Create a script file (myscript) similar to the one below

```
#! /bin/bash
# myscript
set -e

module load cellsnp-lite || exit 1
cd /data/$USER/test/
cellsnp-lite -p $SLURM_CPUS_PER_TASK -s file.bam -O OurDir -R file.csv.gz -b file.tsv --minCOUNT 20
```

2. Submit the script on biowulf: 
 
 
```
[biowulf]$ sbatch --mem=5g -cpus-per-task=4 myscript
```


Using Swarm

Using the 'swarm' utility, one can submit many jobs to the cluster to run concurrently. 


Set up a swarm command file (eg /data/$USER/cmdfile). 



```

cd /data/$USER/dir1; cellsnp-lite ...
cd /data/$USER/dir2; cellsnp-lite ...
cd /data/$USER/dir3; cellsnp-lite ...
...
cd /data/$USER/dir20; cellsnp-lite ...

```


 submit the swarm job:
 
```
$ swarm -f cmdfile --module cellsnp-lite -g 5 -t 4
```

For more information regarding running swarm, see <swarm.html>





















