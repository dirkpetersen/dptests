

document.querySelector('title').textContent = 'DiffIR: Identifying differential intron retention from RNA-seq ';
**DiffIR: Identifying differential intron retention from RNA-seq** 


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



iDiffIR is a tool for identifying differential IR from RNA-seq data. 
It accepts any sorted, indexed BAM file for single- or paired-end reads.



Documentation
* [iDiffIR Github paage](https://github.com/bio-comp/idiffir)
* [iDiffIR Tutorial](https://combi.cs.colostate.edu/idiffir/tutorial.html )


Important Notes
* Module Name: iDiffIR (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Implemented as a Singularity container
* Unusual environment variables set
	+ **IDIFFIR\_HOME**  installation directory
	+ **IDIFFIR\_BIN**    executables directory
	+ **IDIFFIR\_SRC**    source directory
	+ **IDIFFIR\_DATA**  sample data directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive**
[user@cn3316 ~]$ **module load iDiffIR** 
[+] Loading singularity  3.8.5-1  on cn4193
[+] Loading iDiffIR  20220121

```

The iDiffIR application involves a number of executables:
[user@biowulf]$ **ls $IDIFFIR\_BIN** 
build\_classifiers.py getDepths.py realignment\_pipeline.py
classify\_sites.py get\_gene\_expression.py run\_miso\_ir.py
convert\_models.py get\_good\_pairs.py sam\_collate.py
convertSam.py gtf2gff.py sam\_filter.py
ests\_to\_splicegraph.py idiffir\_plotter.py sam\_split.py
f2py idiffir.py sam\_to\_depths.py
find\_splice\_forms.py isolasso\_pipeline.py select\_model\_parameters.py
fix\_unresolved.py isolasso\_update\_graphs.py shell
gene\_model\_to\_splicegraph.py make\_MISO\_AS\_GFF.py simulate\_IR.py
generate\_known\_junctions.py make\_MISO\_IR\_GFF.py smtpd.py
generate\_predicted\_junctions.py plotter.py splicegraph\_statistics.py
generate\_putative\_sequences.py predict\_graphs.py splice\_junction\_pipeline.py
generate\_roc.py predict\_splicegraph.py view\_splicegraph\_multiplot.py
generate\_splice\_site\_data.py psginfer\_pipeline.py
genewise\_statistics.py psginfer\_update\_graphs.py

Their basic usage is as follows: 

```

[user@cn3123 user]$ **idiffir.py -h**
usage: idiffir.py [-h] [-v] [-n] [-l FACTORLABEL FACTORLABEL] [-o OUTDIR] [-s]
                  [-k KRANGE KRANGE] [-c COVERAGE] [-d DEXPTHRESH] [-p PROCS]
                  [-f FDRLEVEL] [-g GRAPHDIRS] [-G GRAPHDIRSONLY]
                  [-m {BF,BH,QV}] [-e {IR,SE}]
                  genemodel factor1bamfiles factor2bamfiles

Identify differentially expressed introns.

positional arguments:
  genemodel             gene model file: NAME.gtf[.gz] | NAME.gff[.gz]
  factor1bamfiles       colon-separated list of bamfiles: PATH-TO-REPLICATE_1
                        [:PATH-TO-REPLICATE_2,...]
  factor2bamfiles       colon-separated list of bamfiles: PATH-TO-REPLICATE_1
                        [:PATH-TO-REPLICATE_2,...]

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose         verbose output [default is quiet running]
  -n, --noplot          Do not plot figures [default is to make figures]
  -l FACTORLABEL FACTORLABEL, --factorlabel FACTORLABEL FACTORLABEL
                        factor labels, example: -f Mutant Wildtype
  -o OUTDIR, --output-dir OUTDIR
                        output file directory name
  -s, --shrink_introns  shrink introns for depth plots [default is no
                        shrinking]
  -k KRANGE KRANGE, --krange KRANGE KRANGE
                        kmin kmax; [default is to search for kmax]
  -c COVERAGE, --coverage COVERAGE
                        coverage cutoff, default = 0.99
  -d DEXPTHRESH, --dexpThresh DEXPTHRESH
                        differential gene expression threshold, [default = 10]
  -p PROCS, --procs PROCS
                        Number of processing cores to use, [default = 1]
  -f FDRLEVEL, --fdrlevel FDRLEVEL
                        FDR test level, [default = 0.05]
  -g GRAPHDIRS, --graph-dirs GRAPHDIRS
                        colon-separated list of directories to recursively
                        search for SpliceGrapher predictions
  -G GRAPHDIRSONLY, --graph-dirs-only GRAPHDIRSONLY
                        colon-separated list of directories to recursively
                        search for SpliceGrapher predictions. In this case
                        only the predicted graphs are used. i.e. the gene
                        models are only used in plots and not for building
                        reduced gene models. Useful for poorly annotated
                        genomes.
  -m {BF,BH,QV}, --multTest {BF,BH,QV}
                        Multiple testing adjustment method BF: Bonferroni, BH:
                        Benjamini-Hochberg, QV: q-values [default = QV]
  -e {IR,SE}, --event {IR,SE}
                        AS event to test, IR: Intron Retention, SE: Exon
                        Skipping [default = IR] [default is IR]

[user@cn3123 user]$ **idiffir\_plotter.py -h**
usage: idiffir_plotter.py [-h] [-v] [-l FACTORLABEL FACTORLABEL] [-o OUTDIR]
                          [-s] [-g GRAPHDIRS] [-p PROCS]
                          genemodel genelist factor1bamfiles factor2bamfiles

Plot highlighted regions between two factors.

positional arguments:
  genemodel             gene model file: NAME.gtf[.gz] | NAME.gff[.gz]
  genelist              File containing gene records to plot. Format for lines
                        is geneID start_1,end_1;...start_n,end_n.
  factor1bamfiles       colon-separated list of bamfiles: PATH-TO-REPLICATE_1
                        [:PATH-TO-REPLICATE_2,...]
  factor2bamfiles       colon-separated list of bamfiles: PATH-TO-REPLICATE_1
                        [:PATH-TO-REPLICATE_2,...]

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose         verbose output [default is quiet running]
  -l FACTORLABEL FACTORLABEL, --factorlabel FACTORLABEL FACTORLABEL
                        factor labels, example: -f Mutant Wildtype
  -o OUTDIR, --output-dir OUTDIR
                        output file directory name
  -s, --shrink_introns  shrink introns for depth plots [default is no
                        shrinking]
  -g GRAPHDIRS, --graph-dirs GRAPHDIRS
                        colon-separated list of directories to recursively
                        search for SpliceGrapher predictions
  -p PROCS, --procs PROCS
                        Number of processing cores to use, [default = 1]

[user@cn3123 user]$ **getDepths.py -h**
usage: getDepths.py [-h] [-o OUTDIR] [-v] bamfile_in

Get chromosomal read depths and junctions

positional arguments:
  bamfile_in            Name of sorted, indexed BAM file

optional arguments:
  -h, --help            show this help message and exit
  -o OUTDIR, --output-dir OUTDIR
                        output file directory name
  -v, --verbose         verbose output

[user@cn3123 user]$ **convertSam.py -h**
usage: convertSam.py [-h] [-o BAMFILE] [-p PROCS] [-m MEMORY] [-v] samfile

Generate sorted BAM and index files for given SAM file

positional arguments:
  samfile               Samfile to convert

optional arguments:
  -h, --help            show this help message and exit
  -o BAMFILE, --outfile BAMFILE
                        Name of converted BAM file [default=.bam]
 -p PROCS, --procs PROCS
 Number of processors to use for BAM sorting (default
 1)
 -m MEMORY, --memory MEMORY
 Max memory (in GBs) for each processor used for BAM
 sorting (default 2)
 -v, --verbose Print verbose output


etc.


```


```

[user@cn3316 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```





