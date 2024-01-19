

document.querySelector('title').textContent = 'Genometools on Biowulf';
Genometools on Biowulf


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



The GenomeTools genome analysis system is a free collection of bioinformatics tools (in the realm of genome informatics) combined into a single binary named gt. It is based on a C library named libgenometools which contains a wide variety of classes for efficient and convenient implementation of sequence and annotation processing software.


Documentation
* <https://github.com/genometools/genometools>


Important Notes
* Module Name: genometools (see [the modules page](/apps/modules.html) for more information)



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive --mem=5g**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load genometools**

[user@cn3144 ~]$ **gt -help**
Usage: bin/gt [option ...] [tool | script] [argument ...]
The GenomeTools genome analysis system.

-i       enter interactive mode after executing 'tool' or 'script'
-q       suppress warnings
-test    perform unit tests and exit
-seed    set seed for random number generator manually.
         0 generates a seed from current time and process id
-help    display help and exit
-version display version information and exit

Tools:

bed_to_gff3
cds
chain2dim
chseqids
clean
...
...

[user@cn3144 ~]$ **gt bed\_to\_gff3 -help**
Usage: bin/gt bed_to_gff3 [BED_file]
Parse BED file and convert it to GFF3.

-featuretype Set type of parsed BED features
             default: BED_feature
-thicktype   Set type of parsed thick BED features
             default: BED_thick_feature
-blocktype   Set type of parsed BED blocks
             default: BED_block
-o           redirect output to specified file
             default: undefined
-gzip        write gzip compressed output file
             default: no
-bzip2       write bzip2 compressed output file
             default: no
-force       force writing to output file
             default: no
-help        display help and exit
-version     display version information and exit

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. gt.sh). For example:



```

#!/bin/bash
set -e
module load genometools
cd /data/$USER
gt bed_to_gff3 -force yes -o out.gff3 input.bed

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --mem=5g gt.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. gt.swarm). For example:



```

cd dir1; gt bed_to_gff3 -force yes -o out.gff3 input.bed
cd dir2; gt bed_to_gff3 -force yes -o out.gff3 input.bed
...
cd dir10; gt bed_to_gff3 -force yes -o out.gff3 input.bed

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f gt.swarm -g 5 --module gt
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module gt Loads the gt module for each subjob in the swarm 
 | |
 | |
 | |












