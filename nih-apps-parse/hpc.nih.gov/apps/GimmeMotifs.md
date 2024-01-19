

document.querySelector('title').textContent = 'GimmeMotifs: de novo motif prediction for ChIP-sequencing experiments ';
**GimmeMotifs: de novo motif prediction for
ChIP-sequencing experiments** 


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



GimmeMotifs is a pipeline for transcription factor motif analysis written in Python.
It incorporates an ensemble of computational tools to predict
motifs de novo from ChIP-sequencing data. Similar
redundant motifs are compared using the weighted information
content similarity score and clustered using an iterative
procedure. A comprehensive output report is generated with several
different evaluation metrics to compare and evaluate the results.



### References:


* Simon J. van Heeringen∗ and Gert Jan C. Veenstra,   
 *GimmeMotifs: a de novo motif prediction pipeline for
ChIP-sequencing experiments.*   

[Bioinformatics](https://academic.oup.com/bioinformatics/article/27/2/270/285575), 2011 **27**(2), 270-271.
* Niklas Bruse, Simon J. van Heeringen  
 *GimmeMotifs: an analysis framework for
transcription factor motif analysis.*   

[bioRxiv, 2018](https://www.biorxiv.org/content/early/2018/11/20/474403) doi: https://doi.org/10.1101/474403


Documentation
	+ [GimmeMotifs documentation](https://gimmemotifs.readthedocs.io/en/master/index.html)
	+ [GimmeMotifs Github page](https://github.com/vanheeringen-lab/gimmemotifs)Important Notes
	+ Module Name: GimmeMotifs (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
	+ Multithreaded
	+ Unusual environment variables set
		- **GMOTIFS\_HOME**  installation directory
		- **GMOTIFS\_BIN**       executable directory
		- **GMOTIFS\_DATA**  sample data directory
Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=16g --cpus-per-task=16**
[user@cn3335 ~]$**module load gimmemotifs** 
[+] Loading GimmeMotifs  0.17.2

```

First, specify genome reference to be used. 
To do this, create a default configuration file

```

[user@cn3335 ~]$**genomepy config generate**
Created config file /home/user/.config/genomepy/genomepy.yaml

```

  

Download sample data to the current folder:

```

[user@cn3335 ~]$**cp -r /usr/local/apps/gimmemotifs/0.17.2/sample\_data/ .**

```

We are now ready to predict new motifs.

```

[user@cn3335 ~]$ **gimme motifs -h**
usage: gimme [-h]  [options] motifs [-h] [-b BACKGROUND]
 [-g GENOME] [--denovo]
 [--known] [--noreport]
 [--rawscore] [--nogc] [-N INT]
 [-p PFMFILE] [-t N]
 [-a ANALYSIS] [-k] [-S]
 [-f FRACTION] [-s N]
 INPUT OUTDIR

positional arguments:
 INPUT FASTA, BED, narrowPeak or region file.
 OUTDIR Output directory.

options:
 -h, --help show this help message and exit
 -b BACKGROUND, --background BACKGROUND
 Background type (random,genomic,gc,promoter,custom) or
 a file with background sequences (FASTA, BED or
 regions)
 -g GENOME Genome name or fasta file
 --denovo Only use de novo motifs
 --known Only use known motifs
 --noreport Don't create a HTML report.
 --rawscore Don't z-score normalize motif scores
 --nogc Don't use GC% bins
 -N INT, --nthreads INT
 Number of threads (default 12)

optional arguments for known motifs:
 -p PFMFILE PFM file with motifs.(default:
 gimme.vertebrate.v5.0.pfm)

optional arguments for de novo motifs:
 -t N, --tools N Tools to use, any combination of AMD,BioProspector,ChI
 PMunk,DiNAMO,GADEM,HMS,Homer,Improbizer,MDmodule,MEME,
 MEMEW,MotifSampler,Posmo,ProSampler,Trawler,Weeder,XXm
 otif (default BioProspector,Homer,MEME)
 -a ANALYSIS, --analysis ANALYSIS
 Analysis type: small, medium, large, xl (xl)
 -k, --keepintermediate
 Don't delete intermediate files
 -S, --singlestrand Only predict motifs for single + strand (default is
 both)
 -f FRACTION, --fraction FRACTION
 Fraction of peaks to use for motif predicton (0.2)
 -s N, --size N Region size to use for motif prediction (200). Set to
 0 to use the size of the input regions.

```

Run 'quick' analysis to predict de novo motifs using the gimme motifs command on the sample data with 16 threads:

```

[user@cn3335 ~]$ **gimme motifs TAp73alpha.bed -g /fdb/genome/hg19 -N 16 -a small -t Homer,MDmodule,BioProspector .** 
2023-01-10 12:41:29,647 - INFO - Creating new config.
2023-01-10 12:41:29,651 - INFO - Using included version of AMD.
2023-01-10 12:41:29,651 - INFO - Using included version of BioProspector.
2023-01-10 12:41:29,651 - INFO - Using included version of ChIPMunk.
2023-01-10 12:41:29,651 - INFO - Using system version of DiNAMO.
...
2023-01-10 12:47:39,662 - INFO - de novo finished
2023-01-10 12:47:39,662 - INFO - output dir: .
2023-01-10 12:47:39,662 - INFO - de novo report: ./gimme.denovo.html
2023-01-10 12:47:59,949 - INFO - creating motif scan tables
2023-01-10 13:26:37,701 - INFO - selecting non-redundant motifs
2023-01-10 13:27:59,200 - INFO - selected 14 non-redundant motifs: ROC AUC 0.952, PR AUC 0.788
2023-01-10 13:27:59,251 - INFO - creating BED files with scan results
2023-01-10 13:28:28,847 - INFO - creating statistics report
...
2023-01-10 13:29:05,250 - INFO - gimme motifs final report: ./gimme.motifs.html


```
