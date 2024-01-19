

document.querySelector('title').textContent = 'Blat on Biowulf';
Blat on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Easyblat](#easyblat)
[Replicating webblat](#webblat)
[Interactive job](#int) 
[Batch job](#sbatch) 
[Swarm of jobs](#swarm) 
 |



BLAT is a DNA/Protein Sequence Analysis program written by Jim Kent at UCSC.
It is designed to quickly find sequences of 95% and greater similarity of
length 40 bases or more. It may miss more divergent or shorter sequence
alignments. It will find perfect sequence matches of 33 bases, and sometimes
find them down to 22 bases. BLAT on proteins finds sequences of 80% and greater
similarity of length 20 amino acids or more. In practice DNA BLAT works well on
primates, and protein blat on land vertebrates.



References & Documentation
* [the BLAT paper](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC187518/)
* [Jim Kent's website](http://www.cse.ucsc.edu/~kent/)


Important Notes
* Module Name: blat (see [the modules page](/apps/modules.html) for more information)
* Blat is singlethreaded. 
* Several genomes are available in /fdb/igenomes, and additional fasta-format databases available in /fdb/fastadb. 
See [here for a list](/refdb).



Easyblat on Biowulf

[Easyblat](/docs/userguide.html#easyblat) is a convenient command-line interface for running Blat on a large number of query sequences.

Put all
your query sequences into a directory, and then type 'easyblat' at the Biowulf
prompt. You will be prompted for all required parameters.

In the example below, blat parameters are set to display only perfect matches (minIdentity=100) . User input in bold.

```

[user@biowulf]$ **easyblat**

EasyBlat: BLAT wrapper for large numbers of sequences on Biowulf
Enter the directory containing your input sequences > **/data/$USER/project/myseqs**
Enter the directory where you want your blat output to go > **/data/$USER/project/out**

The following genomes are available:
0  -  Bos_taurus
1  -  Caenorhabditis_elegans
2  -  Canis_familiaris
3  -  Danio_rerio
4  -  Drosophila_melanogaster
5  -  Gallus_gallus
6  -  Homo_sapiens
7  -  Mus_musculus
8  -  Pan_troglodytes
9  -  Rattus_norvegicus
10  -  Saccharomyces_cerevisiae
11  -  Sus_scrofa
Enter number  from above (e.g. 7) or full path to custom db (e.g /data/$USER/mygenome.fas)
Database to run against > **6**
Versions available for Homo_sapiens:
hg18
hg19
hg38
Enter version of Homo_sapiens genome from above:
Homo_sapiens version:  > **hg38**
Additional blat parameters (e.g. -q type -stepSize=# -minMatch=#) > **-minIdentity=100**

By default, easyblat will allocate 100GB of memory on each node.
Some jobs may require additional memory.
Enter a memory allocation in GB, or leave blank to accept the default >

Local disk scratch space allocation default allocation is 380 GB.
Some jobs may require additional scratch space.
Enter a scratch space allocation in GB, or leave blank to accept the default >


Date                : 26 Feb 2020 17:05:04
Blat version        : 3.5
Query directory     : /data/$USER/project/myseqs
Fasta files         : 10
Output directory    : /data/$USER/project/out
blat db             : /fdb/igenomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa
blat db size        : 3.0GiB
blat options        :  -minIdentity=100
memory allocation   : 100GiB
lscratch allocation : 380GiB
# jobs              : 1
fasta files per job : 11

---------------------- Easyblat - normal mode --------------------------
/usr/local/bin/swarm  --logdir=/data/$USER/project/out/swarm_logs --silent -f /data/$USER/project/out/scripts/swarmfile -t 32 -g 50 --time=8:00:00 --job-name Eblat-26Feb2020-1704 --gres=lscratch:380
Job 49812161 has been submitted.
```


**One-line Easyblat**

Easyblat options can also be provided on the command-line. Note that all options must be provided, else easyblat will prompt for the missing options.

```

[user@biowulf ]$ **easyblat -q /data/$USER/100n -o /data/$USER/out \
 -d /fdb/igenomes/Homo\_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa \
 --blatopts="-minIdentity=100" -m 80 -l 300**

EasyBlat: BLAT wrapper for large numbers of sequences on Biowulf

Date                : 27 Feb 2020 11:29:19
Blat version        : 3.5
Query directory     : /data/$USER/100n
Fasta files         : 100
Output directory    : /data/$USER/out
blat db             : /fdb/igenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa
blat db size        : 2.9GiB
blat options        :  -minIdentity=100
memory allocation   : 80GiB
lscratch allocation : 300GiB
# jobs              : 1
fasta files per job : 101

---------------------- Easyblat - normal mode --------------------------
/usr/local/bin/swarm  --logdir=/data/$USER/out/swarm_logs --silent -f /data/$USER/out/scripts/swarmfile -t 32 -g 80 --time=8:00:00 --job-name Eblat-27Feb2020-1129 --gres=lscratch:300
Job 49859927 has been submitted.

```

**Running against your own database**
You can run against your own database (any fasta format file) by selecting
'other databases', and then entering the full pathname of the database you want
to search. For example:




```

The following genomes are available:
0  -  Bos_taurus
1  -  Caenorhabditis_elegans
2  -  Canis_familiaris
3  -  Danio_rerio
4  -  Drosophila_melanogaster
5  -  Gallus_gallus
6  -  Homo_sapiens
7  -  Mus_musculus
8  -  Pan_troglodytes
9  -  Rattus_norvegicus
10  -  Saccharomyces_cerevisiae
11  -  Sus_scrofa
Enter number  from above (e.g. 7) or full path to custom db (e.g /data/$USER/mygenome.fas)  **/data/$USER/mygenome/genome.fas**


```


**Easyblat options**

Options can be seen with ''easyblat -h'. 


```

NAME
    easyblat - user-friendly script to run BLAT on many large sequences

SYNOPSIS
    easyblat [options]

OPTIONS
    --help, -h
        Show help message

    -q, --querydir=QUERYDIR
        Directory containing fasta files to be used as queries

    -o, --outdir=OUTDIR
        Directory to be used for output and scripts

    -d, --blatdb=blat_DATABASE
        blat database to use

    -m,--mem=memory
        memory to allocate in GB. (100 GB default, typically works well for a
        3 GB genome database)

    -l,--lscratch=lscratch
        local disk to allocate in GB. (380 GB default)

    --blatopts="blat_OPTIONS"
        additional blat options. For example --blatopts="-minIdentity=100"

    -e, --email=EMAIL_ADDRESS
        send mail to EMAIL_ADDRESS at the end of the job

    -n, --dryrun
        Generate all input files and print swarm command, but don't submit.
        Overrides --hold.

    -x, --hold
        Generate all input files, print swarm command, and submit in a held
        state. (Release a hold with 'scontrol release jobid)

DESCRIPTION
    easyblat will submit a swarm of blat jobs that copy the blat database to
    lscratch and allocate an appropriate amount of memory and lscratch.

    It will set up jobs to process sequences from ~100 fasta files in the
    input directory in each job.

    Most options can either be provided on the command line or are prompted
    interactively.

```


Replicating webblat

Users may want to use the same parameter set that the [UCSC web-based blat server](https://genome.ucsc.edu/cgi-bin/hgBlat) uses. As per the 
[blat documentation](https://genome.ucsc.edu/FAQ/FAQblat.html#blat5), the recommended parameter values are: 

```

blat -stepSize=5 -repMatch=2253 -minScore=20 -minIdentity=0 database.2bit query.fa output.psl 

```


For example, to run a large group of sequences against hg19 on Biowulf, use:

```

blat  -stepSize=5  -repMatch=2253  -minScore=20  -minIdentity=0 \
     /fdb/igenomes/Homo_sapiens/UCSC/hg19/hg19.fa    \
     query.fa  query.psl

```


Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load blat**

[user@cn3144 ~]$ **faToNib gi\_22507416.fasta gi\_22507416.nib**

[user@cn3144 ~]$ **blat /fdb/fastadb/hs\_genome.rna.fas gi\_22507416.nib out.psl**
    Loaded 108585020 letters in 42753 sequences
    Searched 1238 bases in 1 sequences 

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. Blat.sh). Since Blat runs spend significant time reading in the database, it is most efficient to 
* allocate as much memory as the size of the database, so that once read, it can be maintained in memory cache
* run a series of blat runs in the same batch job


For example:



```

#!/bin/bash
set -e

module load blat
blat    /fdb/igenomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa   gi_22507416.nib   out1.psl 
blat    /fdb/igenomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa   gi_22507417.nib   out2.psl 
blat    /fdb/igenomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa   gi_22507418.nib   out3.psl 

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command. To determine the best memory allocation, 
check the size of the fasta-format database. e.g.

```

biowulf% ls -lh /fdb/igenomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa
-rwxrwxr-x 1 helixapp helixapp 2.6G Jun 16  2015 genome.fa

```

Therefore, an appropriate memory allocation for this Blat run would be 8 GB.


```
sbatch --mem=8g Blat.sh
```


Note that easyblat, described above, will automatically select the appropriate memory allocation for your jobs. 

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. Blat.swarm). For example:



```

blat    /fdb/igenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa   gi_22507416.nib   out1.psl 
blat    /fdb/igenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa   gi_22507417.nib   out2.psl 
blat    /fdb/igenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa   gi_22507418.nib   out3.psl 

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f Blat.swarm -g 8  --module blat
```

where


|  |  |  |  |
| --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | --module blat Loads the blat module for each subjob in the swarm 
 | |
 | |






























