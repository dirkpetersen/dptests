

document.querySelector('title').textContent = 'nanopolish on Biowulf';
nanopolish on Biowulf


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



nanopolish is a software package for signal-level analysis of Oxford Nanopore sequencing data. Nanopolish can calculate an improved consensus sequence for a draft genome assembly, detect base modifications, call SNPs and indels with respect to a reference genome and more 



Documentation
* [nanopolish Main Site](https://nanopolish.readthedocs.io/en/latest/index.html)


Important Notes
* Module Name: nanopolish (see [the modules page](/apps/modules.html) for more information)
 * Multithreaded app (use -t option; not available in all commands)
 * Environment variables set 
	+ NANOPOLISH\_TESTDATA* Example files in /usr/local/apps/nanopolish/TEST\_DATA



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive --mem 4g**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load nanopolish**

[user@cn3144 ~]$ **cp $NANOPOLISH\_TESTDATA/\* .**

[user@cn3144 ~]$ **tar -xf ecoli\_2kb\_region.tar.gz**

[user@cn3144 ~]$ **cd ecoli\_2kb\_region/**

[user@cn3144 ~]$ **nanopolish index -d fast5\_files/ reads.fasta**
[readdb] indexing fast5_files/
[readdb] num reads: 112, num reads with path to fast5: 112

[user@cn3144 ~]$ **nanopolish**
error: no command provided
usage: nanopolish [command] [options]
  valid commands: 
    --help
    --version
    call-methylation
    eventalign
    extract
    getmodel
    help
    index
    methyltrain
    phase-reads
    polya
    scorereads
    variants
    vcf2fasta
  for help on given command, type nanopolish command --help

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. nanopolish.sh). For example:



```

#!/bin/bash
set -e
module load nanopolish
nanopolish index -d fast5_files/ reads.fasta
```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=2 --mem=4g nanopolish.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. nanopolish.swarm). For example:



```

nanopolish index -d fast5_files1/ reads1.fasta
nanopolish index -d fast5_files2/ reads2.fasta
nanopolish index -d fast5_files3/ reads3.fasta

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f nanopolish.swarm -g 4 -t 2 --module nanopolish
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module nanopolish Loads the nanopolish module for each subjob in the swarm 
 | |
 | |
 | |








