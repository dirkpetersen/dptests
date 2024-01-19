

document.querySelector('title').textContent = 'MAKER: a portable and easily configurable genome annotation pipeline.';
**MAKER: a portable and easily configurable genome annotation pipeline.**


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



MAKER is a portable and easily configurable genome annotation pipeline. 
Its purpose is to allow smaller eukaryotic and prokaryotic genome projects 
to independently annotate their genomes and to create genome databases. 
MAKER identifies repeats, aligns ESTs and proteins to a genome, 
produces ab-initio gene predictions and automatically synthesizes these data 
into gene annotations having evidence-based quality values.




### References:


* Brandi L.Cantarel, Ian Korf, Sofia M.C.Robb, Genis Parra, Eric Ross,
Barry Moore, Carson Holt, Alejandro Sánchez Alvarado, and Mark Yandell   

*MAKER: An easy-to-use annotation pipeline designed for emerging model organism genomes*   

[Genome Research,](https://genome.cshlp.org/content/18/1/188.short)  
2008, v.18: 188–196 doi: 10.1101/gr.6743907.
* Michael S. Campbell Carson Holt Barry Moore and Mark Yandell,   

*Genome Annotation and Curation Using MAKER and MAKER-P*   

[Curr Protoc Bioinformatics,](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4286374) 2014, v.48: 4.11.1–4.11.39; doi: 10.1002/0471250953.bi0411s48.


Documentation
* [MAKER Home page](http://www.yandell-lab.org/software/maker.html)
* [GC\_specific\_MAKER](https://github.com/Childs-Lab/GC_specific_MAKER/blob/master/train_augustus.sh)
* [The MAKER control files explained](http://weatherby.genetics.utah.edu/MAKER/wiki/index.php/The_MAKER_control_files_explained)
* [MAKER Tutorial for WGS Assembly and Annotation](https://weatherby.genetics.utah.edu/MAKER/wiki/index.php/MAKER_Tutorial_for_WGS_Assembly_and_Annotation_Winter_School_2018)
* [Genome Annotation using MAKER](https://gist.github.com/darencard/bb1001ac1532dd4225b030cf0cd61ce2)
* [Genome Annotation and Curation Using MAKER and MAKER-P](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4286374/)


Important Notes
* Module Name: MAKER (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set
	+ **MAKER\_HOME**  installation directory
	+ **MAKER\_BIN**  executable directory
	+ **MAKER\_SRC**  source code directory
	+ **MAKER\_DATA**  sample data directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=20g -n 14**
[user@cn3316 ~]$ **module load maker**
[+] Loading ucsc 390 on cn2350
[+] Loading BLAT  3.5
[+] Loading bedtools  2.29.0
[+] Loading mpich2 3.2.1  for GCC 4.8.5
[+] Loading blast 2.10.0+  ...
[+] Loading repeatmasker  4.1.0  on cn2350
[+] Loading exonerate, version 2.4.0...
[+] Loading perl 5.18.4 on cn2350
[+] Loading maker  2.31.10
[user@cn3316 ~]$ **maker -h**
MAKER version 2.31.10

Usage:

     maker [options] <maker_opts> <maker_bopts> <maker_exe>

Description:

     MAKER is a program that produces gene annotations in GFF3 format using
     evidence such as EST alignments and protein homology. MAKER can be used to
     produce gene annotations for new genomes as well as update annotations
     from existing genome databases.

     The three input arguments are control files that specify how MAKER should
     behave. All options for MAKER should be set in the control files, but a
     few can also be set on the command line. Command line options provide a
     convenient machanism to override commonly altered control file values.
     MAKER will automatically search for the control files in the current
     working directory if they are not specified on the command line.

     Input files listed in the control options files must be in fasta format
     unless otherwise specified. Please see MAKER documentation to learn more
     about control file  configuration.  MAKER will automatically try and
     locate the user control files in the current working directory if these
     arguments are not supplied when initializing MAKER.

     It is important to note that MAKER does not try and recalculated data that
     it has already calculated.  For example, if you run an analysis twice on
     the same dataset you will notice that MAKER does not rerun any of the
     BLAST analyses, but instead uses the blast analyses stored from the
     previous run. To force MAKER to rerun all analyses, use the -f flag.

     MAKER also supports parallelization via MPI on computer clusters. Just
     launch MAKER via mpiexec (i.e. mpiexec -n 40 maker). MPI support must be
     configured during the MAKER installation process for this to work though


Options:

     -genome|g <file>    Overrides the genome file path in the control files

     -RM_off|R           Turns all repeat masking options off.

     -datastore/         Forcably turn on/off MAKER's two deep directory
      nodatastore        structure for output.  Always on by default.

     -old_struct         Use the old directory styles (MAKER 2.26 and lower)

     -base    <string>   Set the base name MAKER uses to save output files.
                         MAKER uses the input genome file name by default.

     -tries|t <integer>  Run contigs up to the specified number of tries.

     -cpus|c  <integer>  Tells how many cpus to use for BLAST analysis.
                         Note: this is for BLAST and not for MPI!

     -force|f            Forces MAKER to delete old files before running again.
                         This will require all blast analyses to be rerun.

     -again|a            recaculate all annotations and output files even if no
                         settings have changed. Does not delete old analyses.

     -quiet|q            Regular quiet. Only a handlful of status messages.

     -qq                 Even more quiet. There are no status messages.

     -dsindex            Quickly generate datastore index file. Note that this
                         will not check if run settings have changed on contigs

     -nolock             Turn off file locks. May be usful on some file systems,
                         but can cause race conditions if running in parallel.

     -TMP                Specify temporary directory to use.

     -CTL                Generate empty control files in the current directory.

     -OPTS               Generates just the maker_opts.ctl file.

     -BOPTS              Generates just the maker_bopts.ctl file.

     -EXE                Generates just the maker_exe.ctl file.

     -MWAS    <option>   Easy way to control mwas_server for web-based GUI

                              options:  STOP
                                        START
                                        RESTART

     -version            Prints the MAKER version.

     -help|?             Prints this usage statement.

```

In order to use maker, copy sample data to your current folder:

```

[user@cn3316 ~]$ **cp $MAKER\_DATA/\*fasta .** 

```

Create control file templates:

```

[user@cn3316 ~]$ **maker -CTL** 
[user@cn3316 ~]$ **ls \*ctl**
maker_bopts.ctl  maker_exe.ctl  maker_opts.ctl

```

Edit these files properly, by adding a missing information to each line before the symbol '#'.   

Alternatively, you can copy the samples of alraady edited files, e.g.:

```

[user@cn3316 ~]$ **cp $MAKER\_HOME/\*ctl .**

```

Now you can run maker by passing it the edited control files as input:

```

[user@cn3316 ~]$ **mpiexec -np 14 maker maker\_opts.ctl maker\_bopts.ctl maker\_exe.ctl** 
STATUS: Parsing control files...
STATUS: Processing and indexing input FASTA files...
STATUS: Setting up database for any GFF3 input...
A data structure will be created for you at:
/gpfs/gsfs7/users/user/MAKER/genome.maker.output/genome_datastore

To access files for individual sequences use the datastore index:
/gpfs/gsfs7/users/user/MAKER/genome.maker.output/genome_master_datastore_index.log

STATUS: Now running MAKER...
examining contents of the fasta file and run log
...
--Next Contig--

Processing run.log file...
MAKER WARNING: The file genome.maker.output/genome_datastore/46/63/chr13//theVoid.chr13/0/chr13.7.te_proteins%2Efa
sta.repeatrunner
did not finish on the last run and must be erased
#---------------------------------------------------------------------
Now starting the contig!!
SeqID: chr13
Length: 114364328
#---------------------------------------------------------------------


examining contents of the fasta file and run log



--Next Contig--

Processing run.log file...
MAKER WARNING: The file genome.maker.output/genome_datastore/63/A4/chr14//theVoid.chr14/0/chr14.0.te_proteins%2Efa
sta.repeatrunner
did not finish on the last run and must be erased
#---------------------------------------------------------------------
Now starting the contig!!
SeqID: chr14
Length: 107043718
#---------------------------------------------------------------------


examining contents of the fasta file and run log



--Next Contig--

Processing run.log file...
MAKER WARNING: The file genome.maker.output/genome_datastore/F7/6D/chr15//theVoid.chr15/0/chr15.1.te_proteins%2Efa
sta.repeatrunner
did not finish on the last run and must be erased
#---------------------------------------------------------------------
Now starting the contig!!
SeqID: chr15
Length: 101991189
#---------------------------------------------------------------------


examining contents of the fasta file and run log

...


```

End the interactive session:

```

[user@cn3316 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```





