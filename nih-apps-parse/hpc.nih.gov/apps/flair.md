

document.querySelector('title').textContent = 'FLAIR: Full-Length Alternative Isoform analysis of RNA ';
**FLAIR: Full-Length Alternative Isoform analysis of RNA** 


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



FLAIR (Full-Length Alternative Isoform analysis of RNA) is a tool for the correction, 
isoform definition, and alternative splicing analysis of noisy reads. 
It is a computational workflow to identify high-confidence transcripts, 
perform differential splicing event analysis, and differential
isoform analysis. 



Documentation
* [FLAIR GitHub page](https://github.com/BrooksLabUCSC/flair)
* [FLAIR Manual](https://flair.readthedocs.io/en/latest/)


### References:


* A.D.Tang, C.M.Soulette, M.J.van Baren, K.Hart, E.Hrabeta-Robinson1, C.J.Wu & A.N.Brooks   

 *Full-length transcript characterization of SF3B1 mutation in chronic lymphocytic leukemia reveals
downregulation of retained introns*   

[Nature Communications volume 11, Article number: 1438 (2020)](https://www.nature.com/articles/s41467-020-15171-6)


Important Notes
* Module Name: FLAIR (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set
	+ **FLAIR\_HOME**  installation directory
	+ **FLAIR\_BIN**       executable directory
	+ **FLAIR\_SRC**       source code directory
	+ **FLAIR\_DATA**  sample data directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --cpus-per-task=4 --mem=4g --gres=lscratch:10**
[user@cn3200 ~]$ **module load FLAIR/1.6.1**
[+] Loading singularity  3.8.5-1  on cn0883
[+] Loading FLAIR  1.6.1
[user@cn3200 ~]$ **flair.py**
usage: python flair.py  --help
modes: align, correct, collapse, quantify, diffExp, diffSplice
Multiple modules can be run when specified using numbers, e.g.:
python flair.py 1234 ...
[user@cn3200 ~]$ **flair.py align -h**
usage: python flair.py align -g genome.fa -r | [options]

flair-align parse options

positional arguments:
 align

optional arguments:
 -h, --help show this help message and exit
 -o O, --output O output file name base (default: flair.aligned)
 -t T, --threads T minimap2 number of threads (4)
 --junction\_bed JUNCTION\_BED
 annotated isoforms/junctions bed file for splice site-
 guided minimap2 genomic alignment
 --pychopper PYCHOPPER
 specify cdna\_classifier.py here to trim reads prior to
 aligning
 -m M, --minimap2 M path to minimap2 if not in $PATH
 --nvrna specify this flag to use native-RNA specific alignment
 parameters for minimap2
 -sam SAM, --samtools SAM
 samtools executable path if not in $PATH
 -c C, --chromsizes C chromosome sizes tab-separated file, used for
 converting sam to genome-browser compatible psl file
 --psl also output sam-converted psl
 --quality QUALITY minimum MAPQ of read alignment to the genome (1)
 -N N retain at most INT secondary alignments from minimap2
 alignment (0)
 --quiet Suppress progress statements from being printed

required named arguments:
 -r R [R ...], --reads R [R ...]
 FastA/FastQ files of raw reads

Either one of the following arguments is required:
 -g G, --genome G FastA of reference genome, can be minimap2 indexed
 --mm\_index MM\_INDEX minimap2 index .mmi file

```

End the interactive session:

```

[user@cn3200 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```





