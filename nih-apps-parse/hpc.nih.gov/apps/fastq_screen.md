

document.querySelector('title').textContent = 'fastq\_screen on Biowulf';
fastq\_screen on Biowulf


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



A tool for multi-genome mapping and quality control. It allows you to screen a library of sequences in FastQ format against a set of sequence databases so you can see if the composition of the library matches with what you expect. 




Features
* It supports multiple threads.
 * MultiQC can parse it's output.



fastq\_screen generates the following graphs demonstrating the proportion of your library was able to map:


![good output graph](/images/good_sequence_screen.png)
In contrast, poor sequencing results will include results from one or more unexpected species. 


![bad output graph](/images/bad_sequence_screen.png)

### References:


* Wingett SW and Andrews S. *FastQ Screen: A tool for multi-genome mapping and quality control* [version 2; referees: 4 approved]. F1000Research 2018, 7:1338; doi:10.12688/f1000research.15931.2; PubMed PMID: 30254741; PubMed Central PMCID:PMC6124377
 [PubMed](https://www.ncbi.nlm.nih.gov/pubmed/30254741) | 
 [PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6124377/) | 
 [Journal](https://f1000research.com/articles/7-1338/v2)


Documentation
* fastq\_screen Main Site:[Main Site](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/)


Important Notes
* Module Name: fastq\_screen (see [the modules page](/apps/modules.html) for more information)
 * Multithreaded
 * This application produces HTML reports. You can use [hpcdrive](/docs/hpcdrive.html) to view these reports on your local workstation.
* Environment variables set 
	+ $FASTQ\_SCREEN\_TEST\_DATA* Example files in $FASTQ\_SCREEN\_TEST\_DATA



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive --cpus-per-task=8**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load fastq\_screen**
[user@cn3144 ~]$ **mkdir /data/$USER/fastq\_screen**
[user@cn3144 ~]$ **cd /data/$USER/fastq\_screen**
[user@cn3144 ~]$ **cp $FASTQ\_SCREEN\_TEST\_DATA/\* .**
[user@cn3144 ~]$ **fastq\_screen --tag fqs\_test\_dataset.fastq.gz**
Using fastq_screen v0.14.0
Reading configuration from '/usr/local/apps/fastq_screen/0.14.0/bin/fastq_screen_v0.14.0/fastq_screen.conf'
Aligner (--aligner) not specified, but Bowtie2 path and index files found: mapping with Bowtie2
Using '/usr/local/apps/bowtie/2-2.4.1/bin/bowtie2' as Bowtie 2 path
Adding database Human
Adding database Mouse
Adding database Rat
Adding database Drosophila
Adding database Worm
Adding database Yeast
Adding database Arabidopsis
Adding database Ecoli
Adding database rRNA
Adding database MT
Adding database PhiX
Adding database Lambda
Adding database Vectors
Adding database Adapters
Using 8 threads for searches
Option --subset set to 0: processing all reads in FASTQ files
Processing fqs_test_dataset.fastq.gz
Not making data subset
Searching fqs_test_dataset.fastq.gz_temp_subset.fastq against Human
Searching fqs_test_dataset.fastq.gz_temp_subset.fastq against Mouse
Searching fqs_test_dataset.fastq.gz_temp_subset.fastq against Rat
Searching fqs_test_dataset.fastq.gz_temp_subset.fastq against Drosophila
Searching fqs_test_dataset.fastq.gz_temp_subset.fastq against Worm
Searching fqs_test_dataset.fastq.gz_temp_subset.fastq against Yeast
Searching fqs_test_dataset.fastq.gz_temp_subset.fastq against Arabidopsis
Searching fqs_test_dataset.fastq.gz_temp_subset.fastq against Ecoli
Searching fqs_test_dataset.fastq.gz_temp_subset.fastq against rRNA
Searching fqs_test_dataset.fastq.gz_temp_subset.fastq against MT
Searching fqs_test_dataset.fastq.gz_temp_subset.fastq against PhiX
Searching fqs_test_dataset.fastq.gz_temp_subset.fastq against Lambda
Searching fqs_test_dataset.fastq.gz_temp_subset.fastq against Vectors
Searching fqs_test_dataset.fastq.gz_temp_subset.fastq against Adapters
Processing complete
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. fastq\_screen.sh). For example:



```

#!/bin/bash
set -e
module load fastq_screen
fastq_screen --tag fqs_test_dataset.fastq.gz

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=8 --mem=2g fastq_screen.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. fastq\_screen.swarm). For example:



```

cd dir1;fastq_screen --tag R1.fq
cd dir2;fastq_screen --tag R2.fq
cd dir3;fastq_screen --tag R3.fq
cd dir4;fastq_screen --tag R4.fq

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f fastq_screen.swarm [-g #] [-t 8] --module fastq_screen
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module fastq\_screen Loads the fastq\_screen module for each subjob in the swarm
 | |
 | |
 | |








