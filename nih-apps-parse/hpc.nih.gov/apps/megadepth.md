

document.querySelector('title').textContent = 'MEGADEPTH: efficient coverage quantification for BigWigs and BAMs.';
**MEGADEPTH: efficient coverage quantification for BigWigs and BAMs.**


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



Megadepth is a fast tool for quantifying alignments and coverage for BigWig and BAM/CRAM input files, using substantially less memory than the next-fastest competitor. Megadepth can summarize coverage within all disjoint intervals of the Gencode V35 gene annotation for more than 19 000 GTExV8 BigWig files in approximately 1 h using 32 threads. Megadepth is available both as a command-line tool and as an R/Bioconductor package providing much faster quantification compared to the rtracklayer package.



### References:


* Wilks, C, Ahmed, O, Baker, DN, Zhang, D, Collado-Torres, L, Langmead, B (2021).   

 *Megadepth: efficient coverage quantification for BigWigs and BAMs.*    

[Bioinformatics,](https://academic.oup.com/bioinformatics/article/37/18/3014/6162880)


Documentation
* [MEGADEPTH Home page](https://github.com/ChristopherWilks/megadepth)


Important Notes
* Module Name: megadepth (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=20g -n 14**
[user@cn3316 ~]$ **module load megadepth**
[+] Loading megadepth  1.2.0
[user@cn3316 ~]$ **megadepth**
megadepth 1.2.0

BAM and BigWig utility.

Usage:
  megadepth  [options]

Options:
 -h --help Show this screen.
 --version Show version.
 --threads # of threads to do: BAM decompression OR compute sums over multiple BigWigs in parallel
 if the 2nd is intended then a TXT file listing the paths to the BigWigs to process in parallel
 should be passed in as the main input file instead of a single BigWig file (EXPERIMENTAL).
 --prefix String to use to prefix all output files.
 --no-auc-stdout Force all AUC(s) to be written to .auc.tsv rather than STDOUT
 --no-annotation-stdout Force summarized annotation regions to be written to .annotation.tsv rather than STDOUT
 --no-coverage-stdout Force covered regions to be written to .coverage.tsv rather than STDOUT
 --keep-order Output annotation coverage in the order chromosomes appear in the BAM/BigWig file
 The default is to output annotation coverage in the order chromosomes appear in the annotation BED file.
 This is only applicable if --annotation is used for either BAM or BigWig input.



```

End the interactive session:

```

[user@cn3316 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```





