

document.querySelector('title').textContent = ' Blast on Biowulf ';

Blast on Biowulf 



|  |
| --- |
| 
Quick Links
[Easyblast](#easyblast)
[Splitting your query seqences](#split)
[Blast databases](#dbs)
[Swarm of jobs](#swarm)
[Single Blast job](#sbatch)
[Centrally-installed Blast databases](/refdb/)
[Documentation](#doc)
 |


The Basic Local Alignment Search Tool (BLAST) finds regions of local similarity between sequences. 
The program compares nucleotide or protein sequences to sequence databases and calculates the statistical significance of matches. 
BLAST can be used to infer functional and evolutionary relationships between sequences as well as help identify members of gene families.

Blast was developed at NCBI, NIH. ([Blast website](http://blast.ncbi.nlm.nih.gov/Blast.cgi))


BLAST on Biowulf is intended for running a large number of
sequence files, such as hundreds or thousands of query sequences, against the
Blast databases. If you have just a few query sequences, you should use Blast
on [the NCBI website](http://blast.ncbi.nlm.nih.gov/Blast.cgi). 

To see all available versions,
use the [modules commands](modules.html) as in the example below. 


```

[user@biowulf ~]$ **module avail blast**

-------------------------- /usr/local/lmod/modulefiles -------------------------
   blast/2.2.26    blast/2.2.30+ (D)

  Where:
   (D):  Default Module
   
[user@biowulf ~]$  **module load blast**

[user@biowulf ~]$ **module list**

Currently Loaded Modules:
  1) blast/2.2.30+
  
[user@biowulf ~]$ **makeblastdb**
USAGE
  makeblastdb [-h] [-help] [-in input_file] [-input_type type]
    -dbtype molecule_type [-title database_title] [-parse_seqids]
    [-hash_index] [-mask_data mask_data_files] [-mask_id mask_algo_ids]
    [-mask_desc mask_algo_descriptions] [-gi_mask]
    [-gi_mask_name gi_based_mask_names] [-out database_name]
    [-max_file_sz number_of_bytes] [-logfile File_Name] [-taxid TaxID]
    [-taxid_map TaxIDMapFile] [-version]

DESCRIPTION
   Application to create BLAST databases, version 2.2.30+

Use '-help' to print detailed descriptions of command line arguments
========================================================================

```


Loading a Blast module will add the blast programs and database-formatting programs to your path.

Easyblast

IMPORTANT: as of Feb 18, 2020, the Blast v5 databases are the default. The Blast
v4 databases are no longer updated by NCBI. To use the old v4 databases with Easyblast, use the -4 flag


Easyblast is an easy interface to Blast on Biowulf. It is a wrapper script that will prompt you 
for all required parameters, set up your jobs appropriately and submit them to the cluster. You will need
to have all your query sequences in multiple files in a single directory (multiple sequences per file is fine). 

Easyblast will run the latest version of Blast+. The version will be printed at the beginning of the Easyblast run,
and will also appear in the summary file in the output directory, and the actual Blast output files. 


Sample session: (user input in bold): 

```

[user@biowulf ~]$ **easyblast**              

EasyBlast: Blast wrapper for large numbers of sequences on biowulf
Enter the directory containing your input sequences > **/data/$USER/blast\_queries**
Enter the directory where you want your Blast output to go > **/data/$USER/blast\_out**
Output directory is not empty. Delete existing files [y|N]? > **y**
Cleaning up output directory... OK

BLAST programs:
    blastn      nucleotide query sequence against nucleotide database
    blastp      protein query sequence against protein database
    blastx      nucleotide query translated in all 6 reading frames
                against a protein database
    tblastn     protein query sequence against a nucleotide database
                translated in all 6 reading frames
    tblastx     6-frame translations of a nucleotide query sequence
                against the 6-frame translations of a nucleotide database
    rpsblast    protein query sequence against protein domain database
    rpstblastn  nucleotide query sequence against protein domain database
	psiblast    protein query sequence against protein database
Which blast program should be run > **blastn**

The following nucleotide databases are available:
    nt - NCBI nonredundant Genbank+EMBL+DDBJ+PDB (no EST, STS, GSS or HTG)
    pdbnt - sequences of nucleotide structures from the Protein Data Bank
    patnt - Patent nucleotide sequences
    GRCh38.p13 - Human Genome GRCh38.p13 (Nov 2019)
    GRCm38.p6 - Mouse Genome GRCm38.p6 (Aug 2019)
Enter short name from above or full path to custom db
Database to run against > **pdbnt**
Additional Blast parameters (e.g. -task megablast -word_size 8 -ungapped) >

By default, easyblast will allocate 10GB of memory on each node.
Some jobs may require additional memory.
Enter a memory allocation in GB, or leave blank to accept the default >

Local disk scratch space allocation default allocation is 380 GB.
Some jobs may require additional scratch space.
Enter a scratch space allocation in GB, or leave blank to accept the default >


Blast version       : 2.10.0+
fasta directory     : /data/$USER/blast/bench/10n
fasta files         : 10
output directory    : /data/$USER/blast/bench/out
blast program       : blastn
blast db            : /fdb/blastdb/pdbnt
blast db file       : /fdb/blastdb/pdbnt.nsq
blast options       :   -num_threads 2
blast db size       : 0.0GiB
memory allocation   : 10GiB
lscratch allocation : 380GiB
# jobs              : 1
fasta files per job : 11

---------------------- Easyblast - normal mode --------------------------
/usr/local/bin/swarm  --logdir=/data/$USER/blast/bench/out/swarm_logs --silent -f /data/$USER/blast/bench/out/scripts/swarmfile -t 32 -g 10 --time=8:00:00 --job-name EBlast-18Feb2020-1025 --gres=lscratch:380
Job 48858512 has been submitted.

```


When the job completes, the output files will appear in your specified output directory.

**Running against your own Blast database**
To run against your own database, enter the db name with full path at the
Database: prompt. For example:



```

       Database to run against: **/data/$USER/blast\_db/my\_db**

```

Database files have suffixes like .nsq, .nin (nucleotide), .psq, .psi (protein) etc.
You should enter the full path and the database name without the suffix. Thus, if your database
files are called my\_db.nsq etc., enter the database name as my\_db. It is expected that your own database is built with the Blast v5 makedb, unless you use the -4 flag to easyblast.

**One-line Easyblast**

Easyblast can also be run with all parameters on the command line. Note that every parameter must be specified, or easyblast will prompt for the missing ones. 


```

[$USER@biowulf easyblast]$ ./easyblast -q /data/$USER/10n -o /data/$USER/out -d /fdb/blastdb/pdbnt -p blastn  --blastopts=" " -m 20 -l 300

EasyBlast: Blast wrapper for large numbers of sequences on biowulf

Important Change 18 Feb 2020: by default, Easyblast uses blast v5 databases
See https://hpc.nih.gov/apps/Blast.html for details

Blast version       : 2.10.0+
fasta directory     : /data/$USER/10n
fasta files         : 10
output directory    : /data/$USER/out
blast program       : blastn
blast db            : /fdb/blastdb/pdbnt
blast db file       : /fdb/blastdb/pdbnt.nsq
blast options       :    -num_threads 2
blast db size       : 0.0GiB
memory allocation   : 20GiB
lscratch allocation : 300GiB
# jobs              : 1
fasta files per job : 11

---------------------- Easyblast - normal mode --------------------------
/usr/local/bin/swarm  --logdir=/data/$USER/out/swarm_logs --silent -f /data/$USER/out/scripts/swarmfile -t 32 -g 20 --time=8:00:00 --job-name EBlast-27Feb2020-1214 --gres=lscratch:300
Job 49870405 has been submitted.

```


**Easyblast options**

Easyblast options can be seen by typing 'easyblast -h'.

```

NAME
    easyblast - user-friendly script to run Blast on many large sequences

SYNOPSIS
    easyblast [options]

OPTIONS
    --help, -h
        Show help message

    -q, --querydir=QUERYDIR
        Directory containing fasta files to be used as queries

    -o, --outdir=OUTDIR
        Directory to be used for output and scripts

    -b, --blastver=BLAST_VERSION
        Version of blast to use. See `module avail` for available versions

    -p, --blastprog=BLAST_PROGRAM
        Blast program to run

    -d, --blastdb=BLAST_DATABASE
        Blast database to use

    -4, --blastdbv4
        Use blast v4 database. These databases are no longer updated as of Oct 2019.

    --blastopts="BLAST_OPTIONS"
        additional blast options. For example --blastopts="-task blastn"

    -e,--email=EMAIL_ADDRESS
        send mail to EMAIL_ADDRESS at the end of the job

    -n, --dryrun
        Generate all input files and print swarm command, but don't submit. Overrides --hold.

    -x, --hold
        Generate all input files, print swarm command, and submit in a held state

DESCRIPTION
    easyblast will submit a swarm of blast jobs that copy the blast database to lscratch and allocate 
    an appropriate amount of memory and lscratch.

    It will set up jobs to process sequences from ~100 fasta files in the input directory in each job.

    Most options can either be provided on the command line or are prompted interactively.

```


Splitting your query sequences

You can put multiple sequences into each of your input sequence files. However,
if you have only 1 sequence file containing a large number of sequences, this
will run on 1 node, with all the Blast runs performed sequentially. You can better
utilize the parallel processing power of Biowulf by dividing your sequences into
a larger number of files. 

 If your query sequences are all in
one file, and you need to split them into multiple sequence files, there are a
couple of utilities available:
* **split\_fasta**: will split a multisequence fasta-format file into a
desired number of files. Usage:

```


 Split_fasta: to split any large uncompressed fasta file 
 Usage: split_fasta [optional parameters] [dir]file.fas
        -n #      number of split files (default=2)
        -o file   root name of output file (default split#)
        -c #      chunks to write out (default 100 entries)
        -d outdir output directory (default = input directory)
        -z        if input file is .Z or .gz compressed


```

Thus, if a file has 100 sequences, and you want to split it into 5
multisequence files, use

```

split_fasta -n 5 sequence_file

```

will produce 5 files, each containing 100/5=20 sequences. The files will be
called split0, split1,...split4.  


```

split_fasta -n 5 -o oligo sequence_file

```

will produce 5 files, each containing 20 sequences. The files will be called
oligo0, oligo1 ..oligo5.



BLAST Databases

Local copies of the sequence databases used by blast can be found in the
directory **/fdb/blastdb**. These data are a (weekly-updated)
mirror of the **<ftp://ftp.ncbi.nlm.nih.gov/blast/db/>** directory
maintained by NCBI. The old v4 databases, which are no longer updated by NCBI, are available in **/fdb/blastdb/v4** for a limited period. 

[Available Blast databases](/refdb)

Running via swarm

Easyblast uses swarm. If you prefer to use swarm directly, you can use easyblast to set up a swarm command file, 
but submit the swarm yourself. 

Use the '-d' flag to easyblast to set up the swarm command file but not submit it. Example:

```

[user@biowulf ~]$ **easyblast -d**

Easyblast running in debug mode. The swarm command file will be generated, but it is left to you to submit it as desired.

EasyBlast: Blast 2.2.30+ for large numbers of sequences
Enter the directory which contains your input sequences: /data/user/blast/bench/10n

[....usual easyblast run....]

---------- Easyblast running in debug mode, no swarm submitted -------------------
Swarm command file created in /home/user/blast_13485.swarm
Recommended swarm submission command
swarm --singleout -f /home/user/blast_13485.swarm --module blast -g 5

```


You can now edit this swarm command file, or submit it requesting more than the default (5G in the example above)
memory. For example

```

swarm --singleout -f /home/user/blast_13485.swarm --module blast -g 10

```



Of course, you can set up a swarm command file yourself and avoid easyblast altogether, along the following lines:


```

# this file is called blastcmd
#
blastn -db /fdb/blastdb/nt -query /data/user/myseqs/seq1.fas -out /data/user/blastout/seq1.out
blastn -db /fdb/blastdb/nt -query /data/user/myseqs/seq2.fas -out /data/user/blastout/seq2.out
blastn -db /fdb/blastdb/nt -query /data/user/myseqs/seq3.fas -out /data/user/blastout/seq3.out
blastn -db /fdb/blastdb/nt -query /data/user/myseqs/seq4.fas -out /data/user/blastout/seq4.out
[...]

```


Determine the size of the database file (see the section on [Blast jobs and memory](#mem)). Let's assume the database is 8.3 GB. Round upwards to 9 GB. 
You should submit the swarm requesting 9 GB with the '-g 9' flag.   

Submit this swarm with

```

swarm -g 9  -f blastcmd --module blast/2.2.26

```


[Blast Database
Update Status](/refdb) -- status of all Blast databases installed on the system.



Single Blast job

You can also submit a single Blast job. Your batch script would be along the following lines:


```

#!/bin/bash

set -e
module load blast
blastx -db /fdb/blastdb/nr -query /data/user/myseqs/seq1.fas -out /data/user/blastout/seq1.out [...other Blast parameters....]

```


If your query sequence file has multiple sequences, it is important to allocate enough memory so that the entire database fits into the allocated memory. If you do not allocate 
enough memory, the Blast program will need to re-read the database multiple times, which will cause huge I/O load and significantly slow down your job. 
In the example above, the database is 'nr'. Based on the calculations [described below](#mem), the nr database is about 75 GB (Aug 2019). So this job should be submitted requesting
about 80 GB.

```

sbatch --mem=80gb  blast.sh

```



Blast jobs and Memory

If you set up your own Blast jobs rather than using Easyblast, note that it is fastest and most efficient 
if the jobs allocate enough memory for the Blast database you're using. If there is insufficient memory,
the database will be read over and over from disk, which is I/O intensive and will slow down your job. It can also
cause huge filesystem problems if you are running a lot of jobs, and you may get contacted by the Biowulf staff. 

Calculate the memory required by totalling the size of the .nsq (for nucleotide databases) or .psq (for 
protein databases) files. e.g.

```

biowulf% du -sh --total /fdb/blastdb/patnt*.nsq
915M    /fdb/blastdb/patnt.00.nsq
858M    /fdb/blastdb/patnt.01.nsq
809M    /fdb/blastdb/patnt.02.nsq
952M    /fdb/blastdb/patnt.03.nsq
958M    /fdb/blastdb/patnt.04.nsq
75M     /fdb/blastdb/patnt.05.nsq
4.5G    total

```

The total of these 6 files is about 4.5 GB. Swarm jobs should therefore be submitted with the '-g 5' flag, or better yet, '-g 6' for safety.

Easyblast will perform these calculations before submitting the job. 

Documentation

* [Blast+ Command Line User Manual](http://www.ncbi.nlm.nih.gov/books/NBK279690/)* [Options for command-line Blast](http://www.ncbi.nlm.nih.gov/books/NBK279675/)* [Extensive Blast documentation](http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs) at the NCBI website.
* Type 'module load blast', then your desired Blast program with '-h'. e.g. 'module load blast; blastn -h'.

























































































