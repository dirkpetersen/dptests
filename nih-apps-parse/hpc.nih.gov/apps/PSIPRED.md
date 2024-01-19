

document.querySelector('title').textContent = 'PSIPRED on Biowulf';
PSIPRED on Biowulf


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



 PSIPRED generates secondary structure predictions using up to four feed-forward neural networks and output from PSI-BLAST.



### References:


* Jones DT.
 [**Protein secondary structure prediction based on position-specific scoring matrices.**](https://www.ncbi.nlm.nih.gov/pubmed/10493868)
*J. Mol. Biol. (1999) 292: 195-202.*


Documentation
* [README](PSIPRED_README.txt)
* [PSIPRED server](http://bioinf.cs.ucl.ac.uk/psipred/)


Important Notes
* Module Name: PSIPRED (see [the modules page](/apps/modules.html) for more information)
* Multithreaded
* Environment variables set 
	+ PSIPRED\_HOME
	+ PSIPRED\_EXAMPLES



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive --cpus-per-task=4**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load PSIPRED**
[user@cn3144 ~]$ **cp $PSIPRED\_EXAMPLES/example.fasta .**
[user@cn3144 ~]$ **runpsipredplus example.fasta**
Running PSI-BLAST with sequence example.fasta ...
Predicting secondary structure...
Pass1 ...
Pass2 ...
Cleaning up ...
Final output files: example.ss2 example.horiz
Finished.
[user@cn3144 ~]$ **ls**
example.fasta  example.horiz  example.ss  example.ss2

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. PSIPRED.sh). For example:



```

#!/bin/bash
module load PSIPRED
runpsipredplus my.fasta

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] PSIPRED.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. PSIPRED.swarm). For example:



```

runpsipredplus file1.fasta > run1.out
runpsipredplus file2.fasta > run2.out
runpsipredplus file3.fasta > run3.out
runpsipredplus file4.fasta > run4.out

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f PSIPRED.swarm [-g #] [-t #] --module PSIPRED
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module PSIPRED Loads the PSIPRED module for each subjob in the swarm 
 | |
 | |
 | |








