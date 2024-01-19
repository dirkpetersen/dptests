

document.querySelector('title').textContent = 'CutRunTools: a flexible pipeline for CUT&RUN processing and footprint analysis ';
**CutRunTools: a flexible pipeline for CUT&RUN processing and footprint
analysis** 


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



CutRunTools is a pipeline for analysis of data produced by the CUT&ampRUN (Cleavage Under Targets and Release Using Nuclease) technology for high-resolution mapping of DNA binding sites.
It is a flexible, general pipeline that facilitates identification of
chromatin-associated protein binding 
and performs genomic footprinting analysis from antibody-targeted
CutRun primary cleavage data. CutRunTools extracts endonuclease cut site information from
sequences of short read fragments and produces single-locus binding estimates, aggregate motif
footprints, and informative visualizations to support the high-resolution mapping capability of
CutRun



### Reference:


* Qian Zhu, Nan Liu, Stuart H. Orkin, Guo-Cheng Yuan   

*CUT&RUNTools: a flexible pipeline for CUT&RUN processing and footprint analysis*   

[bioRxiv](https://www.biorxiv.org/content/10.1101/529081v1.abstract), doi: https://doi.org/10.1101/529081.

Documentation
	+ [CutRunTools Overview and Usage](https://bitbucket.org/qzhudfci/cutruntools/overview)Important Notes
	+ Module Name: CutRunTools (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
	+ Unusual environment variables set
		- **CUTRUNTOOLS\_HOME**  installation directory
		- **CUTRUNTOOLS\_BIN**       executable directory
		- **CUTRUNTOOLS\_SRC**       source code directory
		- **CUTRUNTOOLS\_DATA**  sample data directory
Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=12g -c5 --gres=lscratch:10**
[user@cn3335 ~]$ **module load CutRunTools/20200629** 
[+] Loading bedops  2.4.35
[+] Loading bedtools  2.27.1
[+] Loading bowtie  2-2.3.4.3
[+] Loading macs  2.1.2
[+] Loading meme  4.12.0  on cn2382
[+] Loading HDF5  1.10.4
[+] Loading NetCDF 4.7.4_gcc9.2.0
[+] Loading picard  2.22.2
[+] Loading gcc  7.3.0  ...
[+] Loading GSL 2.4 for GCC 7.2.0 ...
[+] Loading openmpi 3.0.2  for GCC 7.3.0
[+] Loading ImageMagick  7.0.8  on cn2382
[+] Loading pandoc  2.9.2.1  on cn2382
[+] Loading R 3.5.2
[+] Loading samtools 1.9  ...
[+] Loading trimmomatic  0.36  on cn2382
[+] Loading CutRunTools  20200629
[user@cn3335 ~]$**mkdir -p /data/$USER/CutRunTools && cd /data/$USER/CutRunTools**

```

The processing described below involves both interactive and batch submission steps
and follows the outline presented in the CutRunTools Usage documentation (see above).  

First copy the source code and sample data from a system folder to your current folder:

```

[user@cn3335 ~]$ **cp -r $CUTRUNTOOLS\_SRC/\* .**
[user@cn3335 ~]$ **cp -P $CUTRUNTOOLS\_DATA/\* .**

```


Then customize the configuration file config.json for your needs 
by editing/replacing the lines containing the string "user". The current configuration assumes that the input/output directory for your data processing will be a folder "workdir" in your currecvt directory.



  
  


Now edit properly / customize the configuration file config.json and then validate it by running the command:

```

[user@cn3335 ~]$ **./validate.py config.json**

```

If no error messages are produced, run the command:

```

[user@cn3335 ~]$ **./create\_scripts.py config.json**

```

A folder workdir with a number of subfolders and files in it will be created.   

Change your current working directory to that folder:

```

[user@cn3335 ~]$ **cd workdir**

```

The CutRunTools data processing involves four steps.   
  

**STEP 1: Read trimming and alignment.**   

Run either one of the following commands:

```

[user@cn3335 ~]$ **sbatch ./integrated.sh GATA1\_D7\_30min\_chr11\_R1\_001.fastq.gz**

```

or

```

[user@cn3335 ~]$  **./integrated.sh GATA1\_D7\_30min\_chr11\_R1\_001.fastq.gz**

```

Even though the commands specify the \*\_R1\_001.fastq.gz file as the only input, CutRunTools will actually check that both forward and reverse fastq files are present.
  
  

The following new files will be produced in the subfolders of workdir: 

```

aligned.aug10
└──GATA1_D7_30min_chr11_aligned_reads.bam

trimmed
├── GATA1_D7_30min_chr11_1.paired.fastq.gz
├── GATA1_D7_30min_chr11_1.unpaired.fastq.gz
├── GATA1_D7_30min_chr11_2.paired.fastq.gz
└── GATA1_D7_30min_chr11_2.unpaired.fastq.gz

trimmed3
├── GATA1_D7_30min_chr11_1.paired.fastq.gz
└── GATA1_D7_30min_chr11_2.paired.fastq.gz

```

**STEP 2: BAM processing and peak calling.**

```

[user@cn3335 ~]$ **cd aligned.aug10** 

```

Run either one of the following commands:

```

[user@cn3335 ~]$ **sbatch ./integrated.step2.sh GATA1\_D7\_30min\_chr11\_aligned\_reads.bam** 

```

or

```

[user@cn3335 ~]$  **./integrated.step2.sh GATA1\_D7\_30min\_chr11\_aligned\_reads.bam** 

```

A number of output files will be produced, inclusing the peak files \*broadPeak and \*narrowPeak in folders ../macs2.\* and files .stringent.sort.bed in folders ../seacr.\*, as you can see by running the command:

```

[user@cn3335 ~]$  **ls ../macs2.\*/\*Peak ../seacr\*/\*.stringent.sort.bed** 

```

  

**STEP 3: Motif finding.**   

Run the following commands:

```

[user@cn3335 ~]$ **cd ..** 
[user@cn3335 ~]$ **./run\_step3.sh** 

```

This command will submit 12 jobs to the compute cluster. Upon completion of the jobs, 
a number of new files and subfolders inside of the folders macs2.\* and seacr.\* will be produced. 
Alternatively, you can run any of these jobs interactively. For example:

```

[user@cn3335 ~]$ **cd macs2.broad.aug18**
[user@cn3335 ~]$ **./integrate.motif.find.sh GATA1\_D7\_30min\_chr11\_aligned\_reads\_peaks.broadPeak**
[user@cn3335 ~]$ **cd ../macs2.narrow.aug18**
[user@cn3335 ~]$ **./integrate.motif.find.sh GATA1\_D7\_30min\_chr11\_aligned\_reads\_peaks.narrowPeak**
[user@cn3335 ~]$ **cd ../seacr.aug12**
[user@cn3335 ~]$ **./integrate.motif.find.sh GATA1\_D7\_30min\_chr11\_aligned\_reads\_treat.stringent.sort.bed**

```

etc.   
   

**STEP 4: Motif footprinting.**   

Run the following commands:

```

[user@cn3335 ~]$ **cd /data/$USER/CutRunTools/workdir**
[user@cn3335 ~]$ **./run\_step4.sh** 

```

This command will submit 12 jobs to the compute cluster. Upon completion of the jobs, a number of new files and subfolders inside of the folders macs2.\* and seacr.\* will be produced. 
Alternatively, you can run any of these jobs interactively. For example:   


```

[user@cn3335 ~]$ **cd macs2.broad.aug18**
[user@cn3335 ~]$ **./integrate.footprinting.sh GATA1\_D7\_30min\_chr11\_aligned\_reads\_peaks.broadPeak**
[user@cn3335 ~]$ **cd ../macs2.narrow.aug18**
[user@cn3335 ~]$ **./integrate.footprinting.sh GATA1\_D7\_30min\_chr11\_aligned\_reads\_peaks.narrowPeak**
[user@cn3335 ~]$ **cd ../seacr.aug12**
[user@cn3335 ~]$ **./integrate.footprinting.sh GATA1\_D7\_30min\_chr11\_aligned\_reads\_treat.stringent.sort.bed**

```

etc.   
   

End the interactive session:

```

[user@cn3335 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```
