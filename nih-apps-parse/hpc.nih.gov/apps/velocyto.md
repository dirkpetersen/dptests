

document.querySelector('title').textContent = 'velocyto: a library for the analysis of RNA velocity.';
**velocyto: a library for the analysis of RNA velocity.**


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



Velocyto is a library for the analysis of RNA velocity. It includes a command line tool and an analysis pipeline.



### References:


* Gioele La Manno, Ruslan Soldatov, Amit Zeisel, Emelie Braun, Hannah Hochgerner, Viktor Petukhov,
Katja Lidschreiber, Maria E. Kastriti, Peter Lönnerberg, Alessandro Furlan, Jean Fan, Lars E. Borm, Zehua Liu,
David van Bruggen, Jimin Guo, Xiaoling He, Roger Barker, Erik Sundström, Gonçalo Castelo-Branco, Patrick Cramer,
Igor Adameyko, Sten Linnarsson & Peter V. Kharchenko  

*RNA velocity of single cells*    

[Nature](https://www.nature.com/articles/s41586%E2%80%93018%E2%80%930414%E2%80%936) **560**,
494-498 (2018).


Documentation
* [velocyto GitHub page](https://github.com/velocyto-team/velocyto.py)
* [velocito Documentation](http://velocyto.org/velocyto.py/index.html)


Important Notes
* Module Name: velocyto (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set
	+ **VELOCYTO\_HOME**  installation directory
	+ **VELOCYTO\_BIN**       executable directory
	+ **VELOCYTO\_SRC**       source code directory
	+ **VELOCYTO\_DATA**  sample data directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=4g**
[user@cn0882 ~]$**module load velocyto** 
[+] Loading velocyto  0.17  on cn0882
[+] Loading singularity  3.8.2  on cn0882
[user@cn0882 ~]$**velocyto --help** 
Usage: velocyto [OPTIONS] COMMAND [ARGS]...

Options:
  --version  Show the version and exit.
  --help     Show this message and exit.

Commands:
  run            Runs the velocity analysis outputting a loom file
  run10x         Runs the velocity analysis for a Chromium Sample
  run-dropest    Runs the velocity analysis on DropEst preprocessed data
  run-smartseq2  Runs the velocity analysis on SmartSeq2 data (independent bam file per cell)
  tools          helper tools for velocyto

```

Copy sample data to the current folder:

```

[user@cn0882 ~]$ **cp -r $VELOCYTO\_DATA/\* .**

```

Run velocyto on the sample data:

```

[user@cn0882 ~]$ **velocyto run10x --samtools-memory 4000 -m mm10\_rmsk.gtf SAMPLE\_CELLRANGER\_OUTPUT genes.gtf**
2021-09-20 16:30:10,265 - DEBUG - Using logic: Default
...
2021-09-20 16:30:10,266 - DEBUG - Example of barcode: AAACGGGGTAAGTTCC and cell_id: SCAF463_4W_WT_cTEC:AAACGGGGTAAGTTCC-1
...
2021-09-20 16:30:10,385 - WARNING - Not found cell and umi barcode in entry 1 of the bam file
2021-09-20 16:30:10,385 - WARNING - Not found cell and umi barcode in entry 56 of the bam file
2021-09-20 16:30:10,385 - WARNING - Not found cell and umi barcode in entry 57 of the bam file
...
2021-09-20 16:30:16,188 - DEBUG - Parsing Chromosome 1 strand - [line 0]
2021-09-20 16:30:16,627 - DEBUG - Done with 1- [line 54320]
2021-09-20 16:30:16,627 - DEBUG - Assigning indexes to genes
2021-09-20 16:30:16,633 - DEBUG - Seen 846 genes until now
2021-09-20 16:30:16,633 - DEBUG - Parsing Chromosome 1 strand + [line 54321]
2021-09-20 16:30:16,877 - DEBUG - Done with 1+ [line 107106]
2021-09-20 16:30:16,877 - DEBUG - Assigning indexes to genes
2021-09-20 16:30:16,884 - DEBUG - Seen 1727 genes until now
2021-09-20 16:30:16,884 - DEBUG - Parsing Chromosome 10 strand - [line 107107]
2021-09-20 16:30:17,163 - DEBUG - Done with 10- [line 137495]
2021-09-20 16:30:17,163 - DEBUG - Assigning indexes to genes
2021-09-20 16:30:17,167 - DEBUG - Seen 2300 genes until now
2021-09-20 16:30:17,167 - DEBUG - Parsing Chromosome 10 strand + [line 137496]
2021-09-20 16:30:17,317 - DEBUG - Done with 10+ [line 170976]
2021-09-20 16:30:17,317 - DEBUG - Assigning indexes to genes
2021-09-20 16:30:17,320 - DEBUG - Seen 2907 genes until now
...
2021-09-20 16:30:18,831 - DEBUG - Parsing Chromosome 14 strand + [line 399691]
2021-09-20 16:30:18,930 - DEBUG - Done with 14+ [line 422194]
2021-09-20 16:30:18,930 - DEBUG - Assigning indexes to genes
2021-09-20 16:30:18,933 - DEBUG - Seen 7777 genes until now
2021-09-20 16:30:18,933 - DEBUG - Parsing Chromosome 15 strand - [line 422195]
2021-09-20 16:30:19,056 - DEBUG - Done with 15- [line 449456]
...
2021-09-20 16:39:33,489 - DEBUG - Summarizing the results of intron validation.
2021-09-20 16:39:33,741 - DEBUG - Validated 52245 introns (of which unique intervals 21892) out of 582013 total possible introns (considering each possible transcript models).
2021-09-20 16:39:33,741 - DEBUG - Validated 52245 introns (of which unique intervals 21892) out of 582013 total possible introns (considering each possible transcript models).
...
2021-09-20 16:39:33,788 - DEBUG - Read first 0 million reads
 2021-09-20 16:40:52,528 - DEBUG - Read first 10 million reads
2021-09-20 16:41:02,539 - DEBUG - Counting for batch 1, containing 100 cells and 6053119 reads
2021-09-20 16:44:14,661 - DEBUG - 283859 reads not considered because fully enclosed in repeat masked regions
2021-09-20 16:45:46,512 - DEBUG - Read first 20 million reads
2021-09-20 16:46:01,828 - DEBUG - Counting for batch 2, containing 100 cells and 6282651 reads
2021-09-20 16:49:21,786 - DEBUG - 275079 reads not considered because fully enclosed in repeat masked regions
2021-09-20 16:50:48,232 - DEBUG - Read first 30 million reads
2021-09-20 16:51:12,473 - DEBUG - Counting for batch 3, containing 100 cells and 6248290 reads
...

```

End the interactive session:

```

[user@cn0882 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```





