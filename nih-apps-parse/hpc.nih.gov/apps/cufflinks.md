

document.querySelector('title').textContent = 'cufflinks on Biowulf';
cufflinks on Biowulf


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



Cufflinks assembles transcripts, estimates their abundances, and tests for differential expression and regulation
in RNA-Seq samples. It accepts aligned RNA-Seq reads and assembles the alignments into a parsimonious set of transcripts.
Cufflinks then estimates the relative abundances of these transcripts based on how many reads support each one.
Cufflinks is a collaborative effort between the  [Laboratory for Mathematical and Computational Biology](http://bio.math.berkeley.edu/), led by Lior Pachter at UC Berkeley, Steven Salzberg's group
 at the University of Maryland [Center for Bioinformatics and Computational Biology](http://cbcb.umd.edu/),
 and [Barbara Wold's lab](http://woldlab.caltech.edu/) at Caltech.


Cufflinks is provided under the OSI-approved  [Boost License](http://en.wikipedia.org/wiki/Boost_Software_License)


Illumina has provided the RNA-Seq user community with a set of genome sequence indexes (including Bowtie, Bowtie2, and BWA indexes) as
 well as GTF transcript annotation files called iGenomes. These files can be used with TopHat and Cufflinks to quickly
 perform expression analysis and gene discovery. The annotation files are augmented with the tss\_id and p\_id GTF
 attributes that Cufflinks needs to perform differential splicing, CDS output, and promoter user analysis.






Please note that Cufflinks has entered a low maintenance, low support stage as it is now largely superseded by [StringTie](stringtie.html) which provides the same core functionality (i.e. transcript assembly and quantification), in a much more efficient way.



### References:


* **Cufflinks:**
Cole Trapnell, Brian Williams, Geo Pertea, Ali Mortazavi, Gordon Kwan, Jeltje van Baren, Steven Salzberg, Barbara Wold, Lior Pachter.
**Transcript assembly and quantification by RNA-Seq reveals unannotated transcripts and isoform switching during cell differentiation**
*Nature Biotechnology, 2010*
* **Cufflinks -b:**
Adam Roberts, Cole Trapnell, Julie Donaghey, John L. Rinn, Lior Pachter.
**Improving RNA-Seq expression estimates by correcting for fragment bias**
*Genome Biology, 2011*
* **Cufflinks [RABT mode]:**
Adam Roberts, Harold Pimentel, Cole Trapnell, Lior Pachter.
**Identification of novel transcripts in annotated genomes using RNA-Seq**
*Bioinformatics, 2011*
* **Cuffdiff:**
Cole Trapnell, David Hendrickson, Martin Sauvageau, Loyal Goff, John L. Rinn, Lior Pachter
**Differential analysis of gene regulation at transcript resolution with RNA-seq**
*Nature Biotechnology, 2012*


Documentation
* <http://cole-trapnell-lab.github.io/cufflinks/tools/>


Important Notes
* Module Name: cufflinks (see [the modules page](/apps/modules.html) for more information)
* Multithreaded
 * Reference data in /fdb/igenomes/, downloaded from [iGenomes](http://ccb.jhu.edu/software/tophat/igenomes.shtml)


There is a patched version of cufflinks available:



```
module load cufflinks/2.2.1_patched
```

The patch significantly
[accelerates progress at positions where thousands of mate pairs have the same location](https://groups.google.com/forum/#!topic/tuxedo-tools-users/UzLCJhj3lUE). The patched version seems to help when working with the Ensembl human annotation.


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

[user@cn3144 ~]$ **module load cufflinks**
[user@cn3144 ~]$ **cufflinks file.sam**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. cufflinks.sh). For example:



```

#!/bin/bash
cd /data/$USER/mydir
module load cufflinks
cufflinks -p $SLURM_CPUS_PER_TASK inputFile

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] cufflinks.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. cufflinks.swarm). For example:



```

cd /data/$USER/mydir1; cufflinks -p $SLURM_CPUS_PER_TASK inputFile
cd /data/$USER/mydir2; cufflinks -p $SLURM_CPUS_PER_TASK inputFile
cd /data/$USER/mydir3; cufflinks -p $SLURM_CPUS_PER_TASK inputFile
[...]   

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f cufflinks.swarm [-g #] [-t #] --module cufflinks
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module cufflinks Loads the cufflinks module for each subjob in the swarm 
 | |
 | |
 | |








