

document.querySelector('title').textContent = 'mothur on Biowulf';
mothur on Biowulf


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


mothur is a tool for analyzing 16S rRNA gene sequences generated on multiple platforms as part of microbial ecology projects.


### References:


* Schloss PD1, Westcott SL, Ryabin T, Hall JR, Hartmann M, Hollister EB, Lesniewski RA, Oakley BB, Parks DH, Robinson CJ, Sahl JW, Stres B, Thallinger GG, Van Horn DJ, Weber CF. [Introducing mothur: open-source, platform-independent, community-supported software for describing and comparing microbial communities. *Appl Environ Microbiol. (2009)*](https://www.ncbi.nlm.nih.gov/pubmed/19801464)
* Westcott SL, Schloss PD. [De novo clustering methods outperform reference-based methods for assigning 16S rRNA gene sequences to operational taxonomic units. *PeerJ. (2015)*](https://www.ncbi.nlm.nih.gov/pubmed/26664811)


Documentation
* [mothur Main Site](http://www.mothur.org/wiki/Main_Page)
* [mothur MiSeq\_SOP Tutorial](http://www.mothur.org/wiki/MiSeq_SOP)


Important Notes
* Module Name: mothur (see [the modules page](/apps/modules.html) for more information)
* Multithreaded
* Environment variables set 
	+ MOTHUR\_HOME
	+ MOTHUR\_EXAMPLES* Example files in $MOTHUR\_EXAMPLES


 (April 2021) In Mothur version 1.45.2, users attempting classification or alignment with BLAST-based, such as with classify.seqs will likely encounter issues where Mothur is unable to find BLAST binaries. We recommend using kmer classification instead while we work to resolve the problem.




Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --cpus-per-task=8 --mem=8g --gres=lscratch:10**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **cd /lscratch/$SLURM\_JOB\_ID**
[user@cn3144 ~]$ **module load mothur**
[user@cn3144 ~]$ **tar zxf $MOTHUR\_EXAMPLES/MiSeq\_SOP.tgz**
[user@cn3144 ~]$ **cd MiSeq\_SOP**
[user@cn3144 ~]$ **mothur stability.batch**
Linux version

Using ReadLine,Boost,HDF5,GSL
mothur v.1.45.2
Last updated: 4/9/21
by
Patrick D. Schloss

Department of Microbiology & Immunology

University of Michigan
http://www.mothur.org

When using, please cite:
Schloss, P.D., et al., Introducing mothur: Open-source, platform-independent, community-supported software for describing and comparing microbial communities. Appl Environ Microbiol, 2009. 75(23):7537-41.

Distributed under the GNU General Public License

Type 'help()' for information on the commands that are available

For questions and analysis support, please visit our forum at https://forum.mothur.org

Type 'quit()' to exit program

[NOTE]: Setting random seed to 19760620.

Batch Mode

Setting environment variable PROCS to $SLURM_CPUS_PER_TASK

mothur > pcr.seqs(fasta=silva.bacteria.fasta, start=11894, end=25319, keepdots=F, processors=8)

Using 8 processors.

...

mothur > quit()


It took 163 seconds to run 23 commands from stability.batch batch file.

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. mothur.sh). For example:



```
#----- This file is MiSeq_SOP.sh -----#
#!/bin/bash

# Set the environment
module load mothur

# Untar the example data files
tar xzf $MOTHUR_EXAMPLES/MiSeq_SOP.tgz
cd MiSeq_SOP

# Run the batch script
mothur stability.batch

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=8 --mem=4g --job-name=MiSeq_SOP --output=MiSeq_SOP.out MiSeq_SOP.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. mothur.swarm). For example:



```

cd run1; mothur run1.batch
cd run2; mothur run2.batch
cd run3; mothur run3.batch
cd run4; mothur run4.batch

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f mothur.swarm [-g #] [-t #] --module mothur
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module mothur Loads the mothur module for each subjob in the swarm 
 | |
 | |
 | |








