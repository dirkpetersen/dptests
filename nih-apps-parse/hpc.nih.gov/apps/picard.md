

document.querySelector('title').textContent = 'Picard on HPC';
Picard on HPC


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


Picard comprises Java-based command-line utilities that manipulate SAM files, 
 and a Java API (SAM-JDK) for creating new programs that read and write SAM 
 files. Both SAM text format and SAM binary (BAM) format are supported.


Picard changed a bit back to a single jar, so currently picard command lines look like this:

```

   java [java opts] -jar $PICARDJARPATH/picard.jar command [options]
   java -jar $PICARDJARPATH/picard.jar --help
```

### Note 1


MarkDuplicates tends to initiate garbage collection threads. It's suggested that users add `-XX:ParallelGCThreads=5`
 in the picard command and request 6 cpus for the command (1 for picard, 5 for garbage collection)



```

java -Xmx???g -XX:ParallelGCThreads=5 -jar $PICARDJARPATH/picard.jar MarkDuplicates ........

```

###  Note 2


For tools that can generate lots of temperatory files (such as FastqToSam),
or when error message 'slurmstepd: Exceeded job memory limit at some point'
appears, it's suggested to include this flag to the picard command:

```

TMP_DIR=/lscratch/$SLURM_JOBID

```

Then submit with: 

```

sbatch --cpus-per-task=6 --mem=?g --gres=lscratch:200

```

or



```

swarm -f swarmfile -t 6 -g ? --gres=lscratch:200

```

Replace ? above with memory in GB 


Documentation
* <http://broadinstitute.github.io/picard/>



Important Notes
* Module Name: picard (see [the modules 
 page](/apps/modules.html) for more information)
* Multithreaded

 environment variables set 
 + `$PICARDJARPATH` - path to the directory holding the jar 
 file(s)
+ `$PICARD_JARPATH` - path to the directory holding the 
 jar file(s)
+ `$PICARDJAR` - for versions > 1.119, the path to the single 
 picard jar file
+ `$PICARD_JAR` - for versions > 1.119, the path to the 
 single picard jar file





Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --cpus-per-task=5**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load picard**[user@cn3144 ~]$ **java -Xmx4g -XX:ParallelGCThreads=5 -jar $PICARDJARPATH/picard.jar command [options]**
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```




Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. batch.sh). For example:



```

#!/bin/bash
set -e
module load picard
java -Xmx4g -XX:ParallelGCThreads=5 -jar $PICARDJARPATH/picard.jar command [options]
```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=6 batch.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. job.swarm). For example:



```
cd dir1; java -Xmx4g -XX:ParallelGCThreads=5 -jar $PICARDJARPATH/picard.jar command [options]
cd dir2; java -Xmx4g -XX:ParallelGCThreads=5 -jar $PICARDJARPATH/picard.jar command [options]
cd dir3; java -Xmx4g -XX:ParallelGCThreads=5 -jar $PICARDJARPATH/picard.jar command [options]
cd dir4; java -Xmx4g -XX:ParallelGCThreads=5 -jar $PICARDJARPATH/picard.jar command [options]

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f job.swarm -t 6 --module picard
```

where
 

|  |  |
| --- | --- |
| -g *#*  | Number of Gigabytes of memory required for each process (1 line in the swarm command file)
  |
| -t *#* | Number of threads/CPUs required for each process (1 line in the swarm command file).
  |
| --module  | Loads the module for each subjob in the swarm  |














