

document.querySelector('title').textContent = 'RepeatMasker on Biowulf';
RepeatMasker on Biowulf


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



RepeatMasker is a program that screens DNA sequences for interspersed repeats
and low complexity DNA sequences. The output of the program is a detailed
annotation of the repeats that are present in the query sequence as well as a
modified version of the query sequence in which all the annotated repeats have
been masked (default: replaced by Ns). On average, almost 50% of a human
genomic DNA sequence currently will be masked by the program.




Documentation
* [RepeatMasker Main Site](http://www.repeatmasker.org/)
* [RepeatMasker Documentation](http://repeatmasker.org/webrepeatmaskerhelp.html)


Important Notes
* Module Name: repeatmasker (see [the modules page](/apps/modules.html) for more information)
* By default, RepeatMasker will start 2 threads for every CPU in the node resulting in **badly overloaded jobs**. Users must prevent this behavior by setting -pa N where N is equal to the number of allocated CPUs divided by 2. In -pa N, "N" must always be greater than or equal to 2; -pa 1 causes RepeatMasker to run as though the -pa option was not specified and to start 2 threads for each CPU. See the examples below for details. 
* RepeatMasker uses the [Repbase](http://www.girinst.org) libraries. Users should be aware of the [Repbase academic license agreement](repeatmasker/repbase_license.html) before using RepeatMasker on the NIH HPC Systems.
* On Biowulf, RepeatMasker has been configured to have [HMMER](hmmer.html) as the default search engine. [RMBlast](rmblast.html) can be used instead by passing -engine rmblast to the RepeatMasker command.
* **Setting a species/clade:** RepeatMasker in v4.1.1 has switched to the [FamDB format](https://github.com/Dfam-consortium/FamDB) for the Dfam database. You may notice RepeatMasker being more strict with regards to what is acceptable for the -species flag. To check for valid names, you can query the database using famdb.py. See famdb.py --help for usage information and below for an example using our copy of the database:
famdb.py -i /fdb/dfam/Dfam.h5 names mammal



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load repeatmasker**
[user@cn3144 ~]$ **RepeatMasker -pa 2 sequence.fasta *# for any sequence file sequence.fasta***
RepeatMasker version open-4.0.7
Search Engine: HMMER [ 3.1b2 (February 2015) ]
Master RepeatMasker Database: /usr/local/apps/repeatmasker/4.0.7/Libraries/Dfam.hmm ( Complete Database: Dfam_2.0 )


analyzing file sequence.fasta
identifying Simple Repeats in batch 1 of 1
identifying full-length ALUs in batch 1 of 1
identifying full-length interspersed repeats in batch 1 of 1
identifying remaining ALUs in batch 1 of 1
identifying most interspersed repeats in batch 1 of 1
identifying Simple Repeats in batch 1 of 1
processing output: 
cycle 1 
cycle 2 
cycle 3 
cycle 4 
cycle 5 
cycle 6 
cycle 7 
cycle 8 
cycle 9 
cycle 10 
Generating output... 
masking
done
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. repeatmasker.sh). For example:



```

#!/bin/sh
set -e
module load repeatmasker
RepeatMasker -engine rmblast -pa 2 sample.fasta

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] repeatmasker.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. repeatmasker.swarm). For example:



```

RepeatMasker -pa 2 sample1.fasta
RepeatMasker -pa 2 sample2.fasta
RepeatMasker -pa 2 sample3.fasta
RepeatMasker -pa 2 sample4.fasta

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f repeatmasker.swarm [-g #] [-t #] --module repeatmasker
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module repeatmasker Loads the repeatmasker module for each subjob in the swarm 
 | |
 | |
 | |








