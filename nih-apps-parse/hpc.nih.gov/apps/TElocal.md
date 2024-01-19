

document.querySelector('title').textContent = 'TElocal: a tool that utilizes both uniquely and ambiguously mapped reads to quantify transposable element expression at the locus level. ';
**TElocal: a tool that utilizes both uniquely and ambiguously mapped reads to quantify transposable element expression at the locus level.** 


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



TElocal: a tool that utilizes both uniquely and ambiguously mapped reads to quantify transposable element expression at the locus level.



Documentation
* [TElocal Github page](https://github.com/mhammell-laboratory/TElocal)


Important Notes
* Module Name: TElocal (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set
	+ **TElocal\_HOME**  installation directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
  


```

[user@biowulf]$ **sinteractive** 
[user@cn0861 ~]$ **module load TElocal** 
[+] Loading singularity  3.8.5-1  on cn4177
[+] Loading TElocal 1.1.2  ...

```


```

[user@cn0861 ~]$ **TElocal -h** 
usage: TElocal [-h] -b RNAseq.bam --GTF genic-GTF-file --TE TE-annot-file [--stranded option] [--mode TE counting mode]
               [--project name] [--sortByPos] [-i iteration] [--maxL maxL] [--minL minL] [-L fragLength]
               [--verbose [verbose]] [--version]

Measuring TE expression per-locus, per-sample.

optional arguments:
  -h, --help            show this help message and exit
  -b RNAseq.bam, --BAM RNAseq.bam
                        An RNAseq BAM file.
  --GTF genic-GTF-file  GTF or db file for gene annotations
  --TE TE-annot-file    locInd file for transposable element annotations
  --stranded option     Is this a stranded library? (forward, no, or reverse). For "first-strand" cDNA libraries (e.g.
                        TruSeq), choose reverse. For "second-strand" cDNA libraries (e.g. QIAseq stranded), choose
                        forward. DEFAULT: no.
  --mode TE counting mode
                        How to count TE: uniq (unique mappers only), or multi (distribute among all alignments).
                        DEFAULT: multi
  --project name        Name of this project. DEFAULT: TElocal_out
  --sortByPos           Alignment file is sorted by chromosome position.
  -i iteration, --iteration iteration
                        number of iteration to run the optimization. DEFAULT: 100
  --maxL maxL           maximum fragment length. DEFAULT:500
  --minL minL           minimum fragment length. DEFAULT:0
  -L fragLength, --fragmentLength fragLength
                        average fragment length for single end reads. For paired-end, estimated from the input alignment
                        file. DEFAULT: for paired-end, estimate from the input alignment file; for single-end, ignored
                        by default.
  --verbose [verbose]   Set verbose level. 0: only show critical message, 1: show additional warning message, 2: show
                        process information, 3: show debug messages. DEFAULT:2
  --version             show program's version number and exit

Example: TElocal -b RNAseq.bam --GTF gene_annotation.gtf --TE TE_annotation.locInd --mode multi


```

etc.


```

[user@cn0861 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226

```





