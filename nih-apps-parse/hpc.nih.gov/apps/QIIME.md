

document.querySelector('title').textContent = 'QIIME2 on Biowulf';
QIIME2 on Biowulf


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



QIIME 2 is a powerful, extensible, and decentralized microbiome analysis package with a focus on data and analysis transparency. QIIME 2 enables researchers to start an analysis with raw DNA sequence data and finish with publication-quality figures and statistical results.



### References:


* J Gregory Caporaso, Justin Kuczynski, Jesse Stombaugh, Kyle Bittinger, Frederic D Bushman, Elizabeth K Costello, Noah Fierer, Antonio Gonzalez Pena, Julia K Goodrich, Jeffrey I Gordon, Gavin A Huttley, Scott T Kelley, Dan Knights, Jeremy E Koenig, Ruth E Ley, Catherine A Lozupone, Daniel McDonald, Brian D Muegge, Meg Pirrung, Jens Reeder, Joel R Sevinsky, Peter J Turnbaugh, William A Walters, Jeremy Widmann, Tanya Yatsunenko, Jesse Zaneveld and Rob Knight. [QIIME allows analysis of high-throughput community sequencing data.](https://www.ncbi.nlm.nih.gov/pubmed/20383131) *Nature Methods*, 2010


Documentation
* [QIIME2 Main Site](https://qiime2.org/)
* [QIIME2 Documentation](https://docs.qiime2.org)


Important Notes
* Module Name: QIIME (see [the modules page](/apps/modules.html) for more information)
* Some QIIME commands can utilize multiple threads with **--p-n-threads**
* Environment variables:
	+ **$QIIME\_EXAMPLES** -- example files for tutorials* Reference data in **/fdb/QIIME/**
* q2studio is **NOT** available due to system conflicts
* Users connecting to Biowulf using PuTTY with X11-forwarding enabled, may encounter an error with
 matplotlib not outputing .qzv files (using qiime2 demux summarize for example). The workaround
 documented [here](https://forum.qiime2.org/t/error-converting-qza-to-qzv-on-cluster/354) is to
 run the following before running qiime2 commands:
 
```

        mkdir -p $HOME/.config/matplotlib/
        echo "backend: Agg" >> $HOME/.config/matplotlib/matplotlibrc
      
```



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session (user input in bold):



```

[user@biowulf]$ **sinteractive**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load QIIME**
[user@cn3144 ~]$ **cp -R $QIIME\_EXAMPLES/qiime2-moving-pictures-tutorial .**
[user@cn3144 ~]$ **cd qiime2-moving-pictures-tutorial**
[user@cn3144 ~]$ **qiime tools import --type EMPSingleEndSequences \
 --input-path emp-single-end-sequences \
 --output-path emp-single-end-sequences.qza**
[user@cn3144 ~]$ **qiime demux emp-single \
 --i-seqs emp-single-end-sequences.qza \
 --m-barcodes-file sample-metadata.tsv \
 --m-barcodes-column BarcodeSequence \
 --output-dir demux-summarize-out \
 --o-per-sample-sequences demux.qza**
[user@cn3144 ~]$ **qiime demux summarize \
 --i-data demux.qza \
 --o-visualization demux.qzv**
[user@cn3144 ~]$ **qiime dada2 denoise-single \
 --i-demultiplexed-seqs demux.qza \
 --p-trim-left 0 \
 --p-trunc-len 120 \
 --o-representative-sequences rep-seqs-dada2.qza \
 --output-dir denoise-single-out \
 --o-table table-dada2.qza**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. QIIME.sh). For example:



```

#!/bin/bash
module load QIIME
cp -R $QIIME_EXAMPLES/qiime2-atacama-tutorial .
cd qiime2-atacama-tutorial

qiime tools import \
  --type EMPPairedEndSequences \
  --input-path emp-paired-end-sequences \
  --output-path emp-paired-end-sequences.qza

qiime demux emp-paired \
  --verbose \
  --m-barcodes-file sample-metadata.tsv \
  --m-barcodes-category BarcodeSequence \
  --i-seqs emp-paired-end-sequences.qza \
  --o-per-sample-sequences demux \
 --p-rev-comp-mapping-barcodes

qiime demux summarize \
  --verbose \
  --i-data demux.qza \
  --o-visualization demux.qzv

qiime dada2 denoise-paired \
  --verbose \
  --p-n-threads $SLURM_CPUS_PER_TASK \
  --i-demultiplexed-seqs demux.qza \
  --o-table table \
  --o-representative-sequences rep-seqs \
  --p-trim-left-f 13 \
  --p-trim-left-r 13 \
  --p-trunc-len-f 150 \
  --p-trunc-len-r 150

qiime feature-table summarize \
  --verbose \
  --i-table table.qza \
  --o-visualization table.qzv \
  --m-sample-metadata-file sample-metadata.tsv

qiime feature-table tabulate-seqs \
  --verbose \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=8 --mem=20g --gres:lscratch:20 --time=4:00:00 QIIME.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit set of independent commands requiring identical resources.
Create a swarmfile (e.g. QIIME.swarm). For example:



```

qiime dada2 denoise-single \
  --p-n-threads $SLURM_CPUS_PER_TASK \
  --p-trim-left 13 \
  --p-trunc-len 150 \
  --i-demultiplexed-seqs fmt-tutorial-demux-1.qza \
  --o-representative-sequences rep-seqs-1.qza \
  --o-table table-1.qza
qiime dada2 denoise-single \
  --p-n-threads $SLURM_CPUS_PER_TASK \
  --p-trim-left 13 \
  --p-trunc-len 150 \
  --i-demultiplexed-seqs fmt-tutorial-demux-2.qza \
  --o-representative-sequences rep-seqs-2.qza \
  --o-table table-2.qza
qiime dada2 denoise-single \
  --p-n-threads $SLURM_CPUS_PER_TASK \
  --p-trim-left 13 \
  --p-trunc-len 150 \
  --i-demultiplexed-seqs fmt-tutorial-demux-3.qza \
  --o-representative-sequences rep-seqs-3.qza \
  --o-table table-3.qza
qiime dada2 denoise-single \
  --p-n-threads $SLURM_CPUS_PER_TASK \
  --p-trim-left 13 \
  --p-trunc-len 150 \
  --i-demultiplexed-seqs fmt-tutorial-demux-4.qza \
  --o-representative-sequences rep-seqs-4.qza \
  --o-table table-4.qza

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f QIIME.swarm -g 10 -t 8 --module QIIME
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g #  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t #  Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module QIIME  Loads the QIIME module for each subjob in the swarm 
 | |
 | |
 | |








