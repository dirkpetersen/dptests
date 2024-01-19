

document.querySelector('title').textContent = 'ChromHMM on Biowulf';
ChromHMM on Biowulf


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


 ChromHMM can segment genomes into different chromatin states by modeling 
re-occuring combinatorial and spatial pattern of various histone modifications
with a multivariate Hidden Markov Model. The resulting segmentations can
be used to annotate genomes. Bed files for immediate visualization
in genome browsers are generated.


ChromHMM automatically computes state enrichments for functional and
annotation datasets (TSSs, exons, ...) which facilitates the biological
characterization of each state.


On Biowulf, ChromHMM can be run in two ways. The full path to the jarfile can be called
like this:



```

java -mx4000M -Djava.awt.headless=true -jar $CHROMHMM_HOME/ChromHMM.jar [ command ] [ options ]
java -mx4000M -Djava.awt.headless=true -jar $CHROMHMM_JAR [ command ] [ options ]

```

The amount of memory is assigned with -mx*[num]*, where 4000 MB is
allocated. The second option *-Djava.awt.headless=true* is required,
unless an X11 display is available. See [here](/docs/connect.html)
for more information about X11 display.


An easier way is to use the wrapper script **ChromHMM.sh**.
This wrapper script includes an additional option to set the amount of
memory:



```

ChromHMM.sh --memory 8g [ command ] [ options ]

```

By default, ChromHMM uses 4gb of memory. To allocate a different amount of
memory, for example 20gb, include **--memory 20g** on the
commandline.


### References:


* Jason Ernst and Manolis Kellis. *ChromHMM: automating chromatin-state 
 discovery and characterization*. Nature Methods 2012(9): 215.
 [Pubmed](http://www.ncbi.nlm.nih.gov/pubmed/22373907)  | 
 [PMC](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3577932/) | 
 [Journal](https://www.nature.com/articles/nmeth.1906)


Documentation
* [Home page](http://compbio.mit.edu/ChromHMM/)
* [Manual](http://compbio.mit.edu/ChromHMM/ChromHMM_manual.pdf)


Important Notes
* Module Name: ChromHMM (see [the modules page](/apps/modules.html) for more information)
* ChromHMM LernModel can use more than one CPU. Make sure to match the number of cpus 
 requested with the number used with the -p option.
* Module sets the `$CHROMHMM_HOME` and `$CHROMHMM_JAR` environment variables
* Example files in `$CHROMHMM_TEST_DATA`



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$
[user@biowulf]$ **sinteractive --mem=8g --cpus-per-task=4 --gres=lscratch:10**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144]$ **cd /lscratch/$SLURM\_JOB\_ID**
[user@cn3144]$ **module load ChromHMM**
[user@cn3144]$ **cp -R $CHROMHMM\_TEST\_DATA/SAMPLEDATA\_HG18 .**
[user@cn3144]$ # train a model with 8 States; leave some buffer between memory given to JVM and
               # allocated memory
[user@cn3144]$ **ChromHMM.sh --memory 7g LearnModel -p $SLURM\_CPUS\_PER\_TASK SAMPLEDATA\_HG18 out 8 hg18**
...
[user@cn3144]$ **cp -r out /data/user/where/you/want/your/chromhmm/results**

[user@cn3144]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. ChromHMM.sh), which uses the input file 'ChromHMM.in'. For example:



```

#! /bin/bash

function fail() {
    echo "$@" >&2
    exit 1
}
module load ChromHMM/1.20 || fail "could not load module ChromHMM"
cd /lscratch/$SLURM_JOB_ID || fail "could not use lscratch"
cp -R $CHROMHMM_TEST_DATA/SAMPLEDATA_HG18 .
ChromHMM.sh --memory 8g LearnModel -p ${SLURM_CPUS_PER_TASK} \
  SAMPLEDATA_HG18 OUTPUTSAMPLE 10 hg18

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=8 --mem=9g ChromHMM.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. ChromHMM.swarm). For example:



```

ChromHMM.sh --memory 8g LearnModel -p ${SLURM_CPUS_PER_TASK} \
  SAMPLEDATA_HG18 OUTPUTSAMPLE8 8 hg18
ChromHMM.sh --memory 8g LearnModel -p ${SLURM_CPUS_PER_TASK} \
  SAMPLEDATA_HG18 OUTPUTSAMPLE10 10 hg18
ChromHMM.sh --memory 8g LearnModel -p ${SLURM_CPUS_PER_TASK} \
  SAMPLEDATA_HG18 OUTPUTSAMPLE12 12 hg18

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f ChromHMM.swarm -g 9 -t 4 --module ChromHMM/1.20
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g #  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t #  Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module ChromHMM  Loads the ChromHMM module for each subjob in the swarm 
 | |
 | |
 | |








