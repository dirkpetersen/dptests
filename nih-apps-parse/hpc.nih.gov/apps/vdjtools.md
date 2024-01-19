

document.querySelector('title').textContent = 'vdjtools on Biowulf';
vdjtools on Biowulf


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



A comprehensive analysis framework for T-cell and B-cell repertoire sequencing data.



Documentation
* [vdjtools Main Site](https://github.com/mikessh/vdjtools)


Important Notes
* Module Name: vdjtools (see [the modules page](/apps/modules.html) for more information)
* Singlethreaded
 * environment variables set (where ${VER} matches the version loaded) 
	+ $VDJ\_PATH=/usr/local/apps/vdjtools/${VER}
	+ $VDJ\_JAR=/usr/local/apps/vdjtools/${VER}/vdjtools-${VER}.jar* Commands using the circlize R library like PlotFancyVJUsage may fail
 with newer R/circlize versions. See [GitHub Issue 139](https://github.com/mikessh/vdjtools/issues/139) for workaround.



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=4g**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load vdjtools**

[user@cn3144 ~]$ **java -Xmx4g -jar $VDJ\_JAR**
VDJtools V1.1.10

Run as $java -jar vdjtools-1.1.10.jar ROUTINE_NAME arguments

[Basic]
CalcBasicStats
CalcSpectratype
CalcSegmentUsage
PlotFancySpectratype
PlotSpectratypeV
PlotFancyVJUsage

[Diversity]
CalcDiversityStats
RarefactionPlot
PlotQuantileStats

[Overlap]
OverlapPair
CalcPairwiseDistances
ClusterSamples
TestClusters
TrackClonotypes

[Preprocessing]
ApplySampleAsFilter
FilterNonFunctional
FilterByFrequency
DownSample
Decontaminate
FilterBySegment
SelectTop

[Operation]
PoolSamples
JoinSamples
(Enrichment) -> Deprecated

[Annotation]
(ScanDatabase) -> moved to VDJdb since 1.0.5, please visit vdjdb.cdr3.net
CalcCdrAaStats
CalcDegreeStats
Annotate
SegmentsToFamilies

[Util]
FilterMetadata
SplitMetadata
Convert
RInstall

* Run with 'discard_scripts' option prior to ROUTINE_NAME to clean up R scripts upon execution

[user@cn3144 ~]$ **java -Xmx4g -jar $VDJ\_JAR CalcBasicStats -h**
usage: CalcBasicStats [options] [sample1 sample2 sample3 ... if -m is not
                      specified] output_prefix
 -h                         display help message
 -m,--metadata  Metadata file. First and second columns should
 contain file name and sample id. Header is
 mandatory and will be used to assign column
 names for metadata.
 -u,--unweighted Will count each clonotype only once, apart
 from conventional frequency-weighted
 histogram.

[user@cn3144 ~]$ **java -Xmx4g -jar $VDJ\_JAR CalcBasicStats sample1 sample2 sample3**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. vdjtools.sh). For example:



```

#!/bin/bash
module load vdjtools
java -Xmx4g -jar $VDJ_JAR CalcBasicStats sample1 sample2 sample3

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --mem=4g vdjtools.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. vdjtools.swarm). For example:



```

java -Xmx4g -jar $VDJ_JAR CalcBasicStats sample1 
java -Xmx4g -jar $VDJ_JAR CalcBasicStats sample2 
java -Xmx4g -jar $VDJ_JAR CalcBasicStats sample3

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f vdjtools.swarm -g 4 --module vdjtools
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module vdjtools Loads the vdjtools module for each subjob in the swarm 
 | |
 | |
 | |








