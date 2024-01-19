

document.querySelector('title').textContent = 'Strelka on Biowulf';
Strelka on Biowulf


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



Strelka is an analysis package designed to detect somatic SNVs and small indels from the aligned sequencing reads of matched tumor-normal samples. 



### References:


* Saunders et al.Strelka: accurate somatic small-variant calling from sequenced tumor–normal sample pairs. Bioinformatics (2015). [Link](https://academic.oup.com/bioinformatics/article/28/14/1811/218573)


Documentation
* [Strelka Main Site](https://github.com/Illumina/strelka)


Important Notes
* Module Name: strelka (see [the modules page](/apps/modules.html) for more information)
* Multithreaded app
* environment variables set 
	+ $STRELKA\_INSTALL\_DIR* Example files in $STRELKA\_TEST\_DATA



Running strelka is a two step process:


1. Configure the workflow with `configureStrelkaWorkflow.pl`.
 Example configuration files for this step can be found in 
 `$STRELKA_INSTALL_DIR/etc`.
2. Run the makefile generated in step 1


 The current version
still does not support parallelization across more than one node (i.e. the
workflow engine used does not support submitting tasks as jobs with SLURM).


Illumina recommends first running [manta](https://hpc.nih.gov/apps/manta.html) on the samples and
supplying manta's candidate indels as input to strelka with the --indelCandidates
command line option. The example below skips this step.


Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --cpus-per-task=10 --mem=10g**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load strelka/2.7.1**

[user@cn3144 ~]$ **cp -r $STRELKA\_TEST\_DATA/demo . # copy example data**

[user@cn3144 ~]$ **cd demo**

[user@cn3144 ~]$ # configure the workflow with 'configureStrelka${type}Workflow.py'
                 # where type is Starling, Germline, or Somatic

[user@cn3144 ~]$ **configureStrelkaSomaticWorkflow.py \**
              **--normalBam=data/NA12892\_dupmark\_chr20\_region.bam \**
              **--tumorBam=data/NA12891\_dupmark\_chr20\_region.bam \**
              **--referenceFasta=data/chr20\_860k\_only.fa \**
              **--runDir demo\_out**

[user@cn3144 ~]$ **demo\_out/runWorkflow.py -m local -j $SLURM\_CPUS\_PER\_TASK**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. strelka.sh). For example:



```

#!/bin/bash
module load strelka/1.0.15 || exit 1

configureStrelkaSomaticWorkflow.py \
  --normalBam=/path/to/normal.bam \
  --tumorBam=/path/to/tumor.bam \
  --referenceFasta=/path/to/hg19.fa \
  --runDir=/path/to/myAnalysis

# -j N allows make to use up to N parallel processes

demo_out/runWorkflow.py -m local -j $SLURM_CPUS_PER_TASK

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=24 --mem=100g strelka.sh
```

Several things must be changed in this script:


* **/path/to/normal.bam** : Alignment (BAM) file for normal reads
* **/path/to/tumor.bam**  : Alignment (BAM) file for tumor reads
* **/path/to/hg19.fa** : Reference (FASTA) file for the alignment
* **/path/to/myAnalysis** : Output directory


Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. strelka.swarm). For example:



```

demo_out1/runWorkflow.py -m local -j $SLURM_CPUS_PER_TASK
demo_out2/runWorkflow.py -m local -j $SLURM_CPUS_PER_TASK
demo_out3/runWorkflow.py -m local -j $SLURM_CPUS_PER_TASK

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f strelka.swarm -g 100 -t 24 --module strelka
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module strelka Loads the strelka module for each subjob in the swarm 
 | |
 | |
 | |








