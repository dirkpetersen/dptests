

document.querySelector('title').textContent = 'PEPATAC: a modular pipeline for ATAC-seq data processing ';
**PEPATAC: a modular pipeline for ATAC-seq data processing** 


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



PEPATAC is a robust pipeline for Assay for Transposase-Accessible Chromatin using sequencing (ATAC-seq) built on a loosely coupled modular framework.
It may be easily applied to ATAC-seq projects of any size, 
from one-off experiments to large-scale sequencing projects. 
It is optimized on unique features of ATAC-seq data to be fast and accurate 
and provides several unique analytical approaches. 



### References:


* Jason P. Smith, M. Ryan Corces, Jin Xu, Vincent P. Reuter, Howard Y. Chang, and Nathan C. Sheffield,   

*PEPATAC: An optimized pipeline for ATAC-seq data analysis with serial alignments*   

[bioRxiv preprint](https://www.biorxiv.org/content/10.1101/2020.10.21.347054v2.full.pdf)  doi: https://doi.org/10.1101/2020.10.21.3470.


Documentation
* [PEPATAC GitHub page](https://github.com/databio/pepatac)
* [PEPATAC Home page](http://code.databio.org/PEPATAC)
* [PEPATAC\_tutorial](http://pepatac.databio.org/en/latest/tutorial/)


Important Notes
* Module Name: PEPATAC (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set
	+ **PEPATAC\_HOME**  installation directory
	+ **PEPATAC\_BIN**       executable directory
	+ **PEPATAC\_SRC**       source code directory
	+ **PEPATAC\_DATA**  sample data directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --cpus-per-task=16 --mem=32g --gres=lscratch:10**
[user@cn3200 ~]$ **module load PEPATAC/0.10.3** 
[+] Loading pepatac  0.10.3
[+] Loading singularity  3.8.5-1  on cn3089

[user@cn3200 ~]$ **pepatac -h**
usage: pepatac.py [-h] [-R] [-N] [-D] [-F] [-T] [--silent] [--verbosity V] [--logdev] [-C CONFIG_FILE]
                  [-O PARENT_OUTPUT_FOLDER] [-M MEMORY_LIMIT] [-P NUMBER_OF_CORES] [-S SAMPLE_NAME] -I INPUT_FILES
                  [INPUT_FILES ...] [-I2 [INPUT_FILES2 [INPUT_FILES2 ...]]] -G GENOME_ASSEMBLY [-Q SINGLE_OR_PAIRED]
                  [--trimmer {trimmomatic,pyadapt,skewer}] [--aligner {bowtie2,bwa}]
                  [--deduplicator {picard,samblaster,samtools}]
                  [--peak-caller {fseq,fseq2,genrich,hmmratac,homer,macs2}] [-gs GENOME_SIZE]
                  [--peak-type {fixed,variable}] [--extend EXTEND] [--frip-ref-peaks FRIP_REF_PEAKS] [--motif] [--sob]
                  [--no-scale] [--prioritize] [--keep] [--noFIFO] [--lite] [--skipqc]
                  [--prealignment-names PREALIGNMENT_NAMES [PREALIGNMENT_NAMES ...]]
                  [--prealignment-index PREALIGNMENT_INDEX [PREALIGNMENT_INDEX ...]] --genome-index GENOME_INDEX
                  --chrom-sizes CHROM_SIZES [--TSS-name TSS_NAME] [--blacklist BLACKLIST] [--anno-name ANNO_NAME]
                  [--search-file SEARCH_FILE] [-V]

PEPATAC version 0.10.3

optional arguments:
  -h, --help            show this help message and exit
  -R, --recover         Overwrite locks to recover from previous failed run
  -N, --new-start       Overwrite all results to start a fresh run
  -D, --dirty           Don't auto-delete intermediate files
  -F, --force-follow    Always run 'follow' commands
  -T, --testmode        Only print commands, don't run
  --silent              Silence logging. Overrides verbosity.
  --verbosity V         Set logging level (1-5 or logging module level name)
  --logdev              Expand content of logging message format.
  -C CONFIG_FILE, --config CONFIG_FILE
                        Pipeline configuration file (YAML). Relative paths are with respect to the pipeline script.
  -O PARENT_OUTPUT_FOLDER, --output-parent PARENT_OUTPUT_FOLDER
                        Parent output directory of project
  -M MEMORY_LIMIT, --mem MEMORY_LIMIT
                        Memory limit for processes accepting such. Default units are megabytes unless specified using
                        the suffix [K|M|G|T].
  -P NUMBER_OF_CORES, --cores NUMBER_OF_CORES
                        Number of cores for parallelized processes
  -S SAMPLE_NAME, --sample-name SAMPLE_NAME
                        Name for sample to run
  -I2 [INPUT_FILES2 [INPUT_FILES2 ...]], --input2 [INPUT_FILES2 [INPUT_FILES2 ...]]
                        Secondary input files, such as read2
  -Q SINGLE_OR_PAIRED, --single-or-paired SINGLE_OR_PAIRED
                        Single- or paired-end sequencing protocol
  --trimmer {trimmomatic,pyadapt,skewer}
                        Name of read trimming program.
  --aligner {bowtie2,bwa}
                        Name of read aligner.
  --deduplicator {picard,samblaster,samtools}
                        Name of deduplicator program.
  --peak-caller {fseq,fseq2,genrich,hmmratac,homer,macs2}
                        Name of peak caller.
  -gs GENOME_SIZE, --genome-size GENOME_SIZE
                        Effective genome size. It can be 1.0e+9 or 1000000000: e.g. human (2.7e9), mouse (1.87e9), C.
                        elegans (9e7), fruitfly (1.2e8). Default:2.7e9
  --peak-type {fixed,variable}
                        Call variable or fixed width peaks. Fixed width requires MACS2.
  --extend EXTEND       How far to extend fixed width peaks up and downstream.
  --frip-ref-peaks FRIP_REF_PEAKS
                        Path to reference peak set (BED format) for calculating FRiP.
  --motif               Perform motif enrichment analysis.
  --sob                 Use seqOutBias to produce signal tracks, incorporate mappability information, and account for
                        Tn5 bias.
  --no-scale            Do not scale signal tracks: Default is to scale by read count. If using seqOutBias, scales by
                        the expected/observed cut frequency.
  --prioritize          Plot cFRiF/FRiF using mutually exclusive priority ranked features based on the order of
                        feature appearance in the feature annotation asset.
  --keep                Enable this flag to keep prealignment BAM files.
  --noFIFO              Do NOT use named pipes during prealignments.
  --lite                Only keep minimal, essential output to conserve disk space.
  --skipqc              Skip FastQC. Useful for bugs in FastQC that appear with some sequence read files.
  --prealignment-names PREALIGNMENT_NAMES [PREALIGNMENT_NAMES ...]
                        Space-delimited list of prealignment genome names to align to before primary alignment.
  --prealignment-index PREALIGNMENT_INDEX [PREALIGNMENT_INDEX ...]
                        Space-delimited list of prealignment genome name and index files delimited by an equals sign
                        to align to before primary alignment. e.g. rCRSd=/path/to/bowtie2_index/.
  --genome-index GENOME_INDEX
                        Path to primary genome index file. Either a bowtie2 or bwa index.
  --chrom-sizes CHROM_SIZES
                        Path to primary genome chromosome sizes file.
  --TSS-name TSS_NAME   Path to TSS annotation file.
  --blacklist BLACKLIST
                        Path to genomic region blacklist file.
  --anno-name ANNO_NAME
                        Path to reference annotation file (BED format) for calculating FRiF.
  --search-file SEARCH_FILE
                        Required for seqOutBias (--sob). Path to tallymer index search file built with the same read
                        length as the input.
  -V, --version         show program's version number and exit

required named arguments:
  -I INPUT_FILES [INPUT_FILES ...], --input INPUT_FILES [INPUT_FILES ...]
                        One or more primary input files
  -G GENOME_ASSEMBLY, --genome GENOME_ASSEMBLY
                        Identifier for genome assembly

```

Described below are the steps to run PEPATAC on the 
[tutorial example](http://pepatac.databio.org/en/latest/tutorial/).   
  

1. Set up folders:

```

[user@cn3200 ~]$ **mkdir pepatac\_tutorial**
[user@cn3200 ~]$ **export TUTORIAL=$PWD/pepatac\_tutorial**
[user@cn3200 ~]$ **cd $TUTORIAL**
[user@cn3200 ~]$ **mkdir data genomes processed templates tools**
[user@cn3200 ~]$ **cd $TUTORIAL/tools**
[user@cn3200 ~]$ **git clone https://github.com/databio/pepatac.git**
[user@cn3200 ~]$ **cd $TUTORIAL/tools/pepatac**
[user@cn3200 ~]$ **git checkout tags/0.10.3**
Previous HEAD position was 7616783... Merge pull request #199 from databio/dev
HEAD is now at 1348557... Merge pull request #181 from databio/dev

```

2. Initialize refgenie and download assets:

```

[user@cn3200 ~]$ **cd $TUTORIAL/tools**
[user@cn3200 ~]$ **export REFGENIE=$TUTORIAL/refgenie\_config.yaml**
[user@cn3200 ~]$ **refgenie init -c $REFGENIE**
[user@cn3200 ~]$ **refgenie pull hg38/fasta hg38/bowtie2\_index hg38/refgene\_anno hg38/ensembl\_gtf hg38/ensembl\_rb**
[user@cn3200 ~]$ **refgenie build hg38/feat\_annotation**
...
### Arguments passed to pipeline:

* `asset_registry_paths`:  `['hg38/feat_annotation']`
*             `assets`:  `None`
*            `command`:  `build`
*        `config_file`:  `refgenie.yaml`
*             `docker`:  `False`
*              `files`:  `None`
*             `genome`:  `None`
*      `genome_config`:  `None`
* `genome_description`:  `None`
*             `logdev`:  `False`
*          `new_start`:  `False`
*          `outfolder`:  `/data/user/pepatac_tutorial/data`
*             `params`:  `None`
*             `recipe`:  `None`
*            `recover`:  `False`
*       `requirements`:  `False`
*             `silent`:  `False`
*     `skip_read_lock`:  `False`
*    `tag_description`:  `None`
*          `verbosity`:  `None`
*            `volumes`:  `None`
...
### Pipeline completed. Epilogue
*        Elapsed time (this run):  0:02:19
*  Total elapsed time (all runs):  0:19:49
*         Peak memory (this run):  0.0739 GB
*        Pipeline completed time: 2022-07-20 10:18:01
Finished building 'feat_annotation' asset
...
[user@cn3200 ~]$ **refgenie pull rCRSd/fasta**
[user@cn3200 ~]$ **refgenie pull rCRSd/bowtie2\_index**

```

Add the export REFGENIE line to your .bashrc or .profile to ensure it persists.   


3. Download tutorial read files:

```

[user@cn3200 ~]$ **wget http://big.databio.org/pepatac/tutorial1\_r1.fastq.gz**
[user@cn3200 ~]$ **wget http://big.databio.org/pepatac/tutorial1\_r2.fastq.gz**
[user@cn3200 ~]$ **wget http://big.databio.org/pepatac/tutorial2\_r1.fastq.gz**
[user@cn3200 ~]$ **wget http://big.databio.org/pepatac/tutorial2\_r2.fastq.gz**

[user@cn3200 ~]$ **mv \*.fastq.gz $TUTORIAL/tools/pepatac/examples/data/**

```

4. Configure project files:

```

[user@cn3200 ~]$ **cd $TUTORIAL/tools/pepatac/examples/tutorial**

```

Edit the files tutorial.csv and tutorial.yaml if needed.   


```

[user@cn3200 ~]$ **cd $TUTORIAL**
[user@cn3200 ~]$ **cp $PEPATAC\_CONFIG/compute\_config.yaml .**
[user@cn3200 ~]$ **export DIVCFG=$TUTORIAL/compute\_config.yaml**

```

Add the export DIVCFG line to your ~/.bashrc to ensure it persists.

```

[user@cn3200 ~]$ **cd $TUTORIAL/templates**
[user@cn3200 ~]$ **cp $PEPATAC\_CONFIG/localhost\_template.sub .**

```

The PEPATAC pipeline is divided into two major parts:   

1) first, it processes each sample individually at the sample level;   

2) once sample processing is complete, the project-level part aggregates, analyzes, and summarizes
the results across samples.   
   

5. Run the PEPATAC pipeline at the **sample level** using looper:

```

[user@cn3200 ~]$ **cd $TUTORIAL/tools/pepatac**
[user@cn3200 ~]$ **looper run -i examples/tutorial/tutorial\_refgenie.yaml**
...
Looper version: 1.3.1
Command: run
...

Looper finished
Samples valid for job generation: 2 of 2
Commands submitted: 2 of 2
Jobs submitted: 2

```

The previous command produces two sbatch scripts.   
   

Finally, submit these scripts to the cluster:

```

[user@cn3200 ~]$  **sbatch $TUTORIAL/processed/submission/PEPATAC\_tutorial1.sub**
[user@cn3200 ~]$  **sbatch $TUTORIAL/processed/submission/PEPATAC\_tutorial2.sub**

```

6. Run the PEPATAC pipeline at the **project level** using looper:

```

[user@cn3200 ~]$ **cd $TUTORIAL/tools/pepatac**
[user@cn3200 ~]$ **looper runp examples/tutorial/tutorial\_refgenie.yaml** 
Looper version: 1.3.2
Command: runp
...
............................................................
............................................................
............................................................
...

```

End the interactive session:

```

[user@cn3200 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```





