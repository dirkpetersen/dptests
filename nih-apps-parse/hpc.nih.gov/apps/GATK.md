

document.querySelector('title').textContent = 'GATK on Biowulf';

GATK on Biowulf

Description
The Genome Analysis Toolkit (GATK) is a software package developed at the
Broad Institute to analyze high-throughput sequencing data. The toolkit
includes a wide variety of tools, with a focus on variant discovery and
genotyping as well as emphasis on data quality assurance.




---


[### Online Tutorial: A practical introduction to GATK 4 on Biowulf](https://hpc.nih.gov/training/gatk_tutorial)
**New in May 2021:** A [self-paced, online tutorial](https://hpc.nih.gov/training/gatk_tutorial) to work through 
a GATK example on Biowulf. 
Developed by the Biowulf staff, this tutorial includes a case study of germline variant discovery with WGS data from a trio, and benchmarks for each step. By working through the tutorial, you will learn NGS data preprocessing and how to optimize your Biowulf batch jobs. 




---


### Notes


* With the release of GATK 3.5-0 MuTect2 and ContEst are
 included as part of GATK.
* GATK 4 is significantly different from GATK 3. Both are documented
 separately below.
* GATK 4 is the default module as of April 2021.
* In most cases, GATK jobs can run successfully on norm partition, there is no need/benefit to run GATK on largemem partition.


There are multiple versions of GATK available. An easy way of selecting the
version is to use [modules](/apps/modules.html). To see the modules
available, type



```

    module avail GATK 

```

To select a module use



```

    module load GATK/[version]

```

where `[version]` is the version of choice.


### Documentation


* [Home page](https://www.broadinstitute.org/gatk/)
* [Guide](https://www.broadinstitute.org/gatk/guide/)
* [Tool Documentation](https://www.broadinstitute.org/gatk/guide/tooldocs/)
* [Forum](http://gatkforums.broadinstitute.org/)
* [Blog](https://www.broadinstitute.org/gatk/blog)


### Release notes


* [3.3 Release notes](http://gatkforums.broadinstitute.org/discussion/4739/release-notes-for-gatk-version-3-3)
* [3.4 Release notes](http://gatkforums.broadinstitute.org/discussion/5562/release-notes-for-gatk-version-3-4)
* [3.5 Release notes](http://gatkforums.broadinstitute.org/discussion/6494/release-notes-for-gatk-version-3-5)
* [3.6 Release notes](http://gatkforums.broadinstitute.org/wdl/discussion/7711/release-notes-for-gatk-version-3-6)
* [3.7 Release notes](http://gatkforums.broadinstitute.org/gatk/discussion/8691/release-notes-for-gatk-version-3-7)
* [3.8 Release notes](https://gatkforums.broadinstitute.org/gatk/discussion/10062/release-notes-for-gatk-version-3-8)
* Release notes for GATK 4 can be found on the [GitHub page](https://github.com/broadinstitute/gatk/releases)


### References


* Aaron McKenna et al. *The Genome Analysis Toolkit: a MapReduce framework for
 analyzing next-generation DNA sequencing data.* Genome Research 2010,
 20:1297-1303.
 [Pubmed](http://www.ncbi.nlm.nih.gov/pubmed/20644199) | 
 [PMC](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2928508/) | 
 [Journal](http://genome.cshlp.org/content/20/9/1297)
* Mark A. DePristo et al. *A framework for variation discovery and genotyping
 using next-generation DNA sequencing data.* Nature Genetics 2011, 43:491-498.
 [Pubmed](http://www.ncbi.nlm.nih.gov/pubmed/21478889) | 
 [PMC](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3083463/) | 
 [Journal](https://www.nature.com/articles/ng.806)
* Geraldine A. Van der Auwera et al. *From FastQ Data to High-Confidence
 Variant Calls: The Genome Analysis Toolkit Best Practices Pipeline.*
 Current Protocols In Bioinformatics 2013, 11:11.10:11.10.1–11.10.33.
 [Pubmed](http://www.ncbi.nlm.nih.gov/pubmed/25431634) | 
 [PMC](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4243306/) | 
 [Journal](http://onlinelibrary.wiley.com/doi/10.1002/0471250953.bi1110s43/abstract;jsessionid=2E0A8F753564C2126232257E21B17766.f04t01)
* Poplin et al. *Scaling accurate genetic variant discovery to tens of thousands of samples*
 bioRxiv,2017
 [bioRxiv](https://doi.org/10.1101/201178) |



Resource bundles
GATK
[Resource](http://gatkforums.broadinstitute.org/discussion/1213/whats-in-the-resource-bundle-and-how-can-i-get-it)
[bundles](https://www.broadinstitute.org/gatk/guide/article.php?id=1213)
contain reference genome assembles (.fa, .fai, .dict), variation data (dbSNP,
HapMap, 1000 genomes, ...), and some genotype calls for a gold standard
genome.


Resource bundles are available at


**`/fdb/GATK_resource_bundle`**



* GATK 3
* GATK 4








|  |
| --- |
| 
Quick Links
[Overview](#overview3)
[Batch job on Biowulf](#serial3)
[Swarm of jobs](#swarm3)
[Interactive job on Biowulf](#int3)
 |




Overview
### Environment module


The environment module for GATK (`module load GATK` or
 `module load GATK/version`)


* loads R as a dependency
* sets the environment variables
	+ `$GATK_HOME`
	+ `$GATK_JAR`
	+ `$GATK_JARPATH`
	+ `$GATK_KEY`
	+ `$GATK_QUEUE_JAR`
	+ `$GATK_TEST_DATA`
* adds a GATK wrapper script developed in house and utilities for lifting
 over VFC files between genome builds


### Wrapper


The GATK wrapper script provides a simplified interface to the GATK tools:



```

usage: GATK [options] [-m MEM] tool [gatk_tool_args]

This GATK wrapper script invokes GATK tools with a more friendly interface.
For documentation on each of the tools and their arguments see
https://www.broadinstitute.org/gatk/guide/tooldocs/

positional arguments:
  tool               GATK tool. Use 'LIST' or 'list' to get the current list
                     of available tools

optional arguments:
  -h, --help         show this help message and exit
  -m MEM, --mem MEM  maximal heap memory allowed for the Java process.
                     *Format*: <number>[kmg]. *Examples*: 1024k, 1000m, 2g,
                     32g. Passed to Java with -Xmx (default: 2g)
  -p N, --pgct N     maximal number of parallel GC threads allowed for Java
                     VM. Defaults to 1 less than the number of allocated CPUs.
                     (default: -1)
  -n, --dry-run      show command that would have been executed and set up
                     temp dir (default: False)

```

It takes care of finding the GATK jar and writing temporary files to either


* `/lscratch/${SLURM_JOB_ID}`: if running on a compute node with
 allocated local disk space
* `./unique_tmp_dir`: if running on a compute 
 node without allocated local disk space


Whenever possible, please use lscratch for temporary files. It will be
 more performant and reduce strain on the shared filesystems.


### Calling GATK directly


When not using the `GATK` wrapper script, the basic GATK command
 should be 



```

java -Xmx####g -Djava.io.tmpdir=/lscratch/${SLURM_JOB_ID} \
  -XX:ParallelGCThreads=## -jar $GATK_JAR -T <Toolname> [options]

```

which requires requesting lscratch space at job submission.


Note that `/lscratch` is the node-local scratch space. See the
 [user guide](/docs/b2-userguide.html#local) for more information
 about allocating local scratch space.


Some tools may perform better with non-default garbage collectors.


#### Memory


The `-Xmx` option to Java specifies the maximal size of the
 memory allocation pool. It accepts suffixes `k`, `m`, and
 `g`. This is the option used by the wrapper script to limit memory
 usage.


Note that you should request additional memory (1-2GB) for batch/interactive jobs on
 top of the memory given to java to ensure that the whole java process stays
 below the allocated memory limit.


#### Redirect temporary files to /lscratch/${SLURM\_JOB\_ID}


GATK writes some temporary files during its run. By default, Java writes
 temp files to `/tmp`. However this is not advisable because
 `/tmp` is used by the operating system and is very limited in space.
 The tmp files should be redirected to `/lscratch/${SLURM_JOB_ID}` (if
 allocating local disk space) 


```

-Djava.io.tmpdir=/lscratch/${SLURM_JOB_ID}

```

to the java command as shown above. The wrapper script takes care of this
 automatically. Note that the temp dir must already exist if
 specifying it manually. The wrapper script will automatically create it
 if necessary.


#### Limit the number of parallel garbage collection threads


The Java VM has a tendency to start many threads including parallel threads
 for garbage collection since it automatically detects all CPUs on a compute
 node irrespective of the allocation. This can be aggravated when there is not
 sufficient memory allocated to the JVM. Use -XX:ParallelGCThreads=##
 to limit the number of parallel garbage collection threads to the number of
 allocated CPUs up to 8 CPUs plus 5 for every 8 beyond that. The GATK wrapper
 by default sets the number of GC threads correctly one less than the number of
 allocated CPUs.


GATK 3.8-0 made the intel de/infaltor the default. While
 slightly faster than their jdk equivalents, this did introduce a bug that can
 lead to segfaults when allocating large heaps. Adding `--use_jdk_inflater
 --use_jdk_deflater` to the GATK tool options avoids this. The GATK
 wrapper takes care of this automatically.


### Phone home


Note that GATK has a
 [phone
 home](https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_engine_CommandLineGATK.php#--phone_home) feature which will send some generic information about each GATK run
 to the developers at the Broad Institute. This feature will not work on the
 Biowulf computational nodes, which are on a private network. This will not
 affect your runs in most cases. If you're getting an error related to
 'phone-home', add  `-et NO_ET -K $GATK_KEY`
 to your GATK command line to disable the phone-home feature. 


### Open files


GATK opens many files simultaneously. It can happen that you hit the
 open-file limit on Biowulf, and see errors like:



```

##### ERROR MESSAGE: Couldn't read file .... (Too many open files)

```

The most important change you can make is to reduce the multithreading. If
 your GATK command has any multithreading flags turned on (e.g. '-nt' and
 'nct'), reduce the number of threads to 1.


If that happens, the limit on the number of open files has to be
 increased from its default soft limit of 4096. This can be done by adding
 the following command to your batch script: `ulimit -n 16384`.


Here are some other possibilities (thanks to Ryan Neff, NHGRI):


* Reduce the number of open files. When doing variant calling, first create gVCF
 files for each sample and then combine them together in batches of 100-200, and
 then do the final variant calls on the remaining 5-20 files (GATK's
 HaplotypeCaller and GenotypeGVCFs). If you are using the UG or a different
 variant caller, try first discovering variant calls on a per-sample basis (-gt
 mode DISCOVERY), combining all of those files together, and then genotyping
 each sample on it's own from the list of all variants (-gt mode
 GENOTYPE\_GIVEN\_ALLELES).
* Reduce the number of open file pointers. This can be done by increasing RAM,
 reducing the number of temporary files open at one time. Sometimes, in the
 GATK, setting --disableAutoIndexCreationAndLockingWhenReadingRods helps with
 performance and other I/O issues.


Relevant discussions from the GATK mailing list:


* [GATK\_UnifiedGenotyper\_Unable to merge temporary Tribble output file](http://gatkforums.broadinstitute.org/discussion/1404/gatk-unifiedgenotyper-unable-to-merge-temporary-tribble-output-file)
* [GATK Unified Genotyper Too many files open even with ulimit](http://gatkforums.broadinstitute.org/discussion/3627/gatk-unified-genotyper-too-many-files-open-even-with-ulimit)


### GATK Parallelization


Some GATK tools can make use of more than one CPU using option `-nt` 
 or `-nct`. See the GATK documentation for more detail.




Running a single GATK batch job on Biowulf
Set up a batch script along the following lines:



```

#! /bin/bash
# gatk.bat
set -e

module load GATK/3.8-0 || exit 1

cd /data/$USER/test_data/GATK

# use 8 GB memory
java -Xmx8g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID -jar $GATK_JAR \
  -T CountReads \
  -R exampleFASTA.fasta \
  -I exampleBAM.bam \
  -et NO_ET -K $GATK_KEY

```

Submit to the queue with [sbatch](/docs/userguide.html):



```

biowulf$ **sbatch --mem=9g --gres=lscratch:10 gatk.bat**

```

Since the simple `CountReads` tool does not support
 multithreading, only a single core (the default) is requested.




Running a swarm of GATK batch jobs on Biowulf
The following swarm commend file would run the [RealignerTargetCreator](https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_indels_RealignerTargetCreator.php)
 ,which determines which genomic intervals require realignment, for a number of
 samples:



```

# this file is gatk-swarm
GATK -m 7g RealignerTargetCreator \
  -R ref.fasta
  -I input1.bam \
  -o output1.intervals \
  --known /fdb/GATK_resource_bundle/hg19/1000G_phase1.indels.hg19.vcf.gz
GATK -m 7g RealignerTargetCreator \
  -R ref.fasta
  -I input2.bam \
  -o output2.intervals \
  --known /fdb/GATK_resource_bundle/hg19/1000G_phase1.indels.hg19.vcf.gz
[...]

```

These command lines specify that the process will use 7 GB of memory
 (-Xmx7g). Thus, this swarm job should be submitted with the '-g 7' flag, to
 tell swarm that each process will require 7 GB of memory.



```

biowulf$ **swarm -g 7 -f gatk-swarm --module GATK/3.8-1**

```



Interactive job on Biowulf
Allocate an interactive session with [sinteractive](/docs/userguide.html#int)
 and use as shown below



```

[user@biowulf]$ **sinteractive --gres=lscratch:50 --cpus-per-task=2 --mem=6g**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144]$ **cd /lscratch/$SLURM\_JOB\_ID**
[user@cn3144]$ **module load GATK/3.8-1**    
[user@cn3144]$ **cp ${GATK\_TEST\_DATA:-none}/NA12878.ba[im]\* .**
[user@cn3144]$ **ls -lh**
total 18M
-rw-r--r-- 1 user group 2.2M Jun 30 13:45 NA12878.bai
-rw-r--r-- 1 user group  15M Jun 30 13:45 NA12878.bam
[user@cn3144]$ **GATK CountReads \
 -R /fdb/GATK\_resource\_bundle/hg38/Homo\_sapiens\_assembly38.fasta \
 -I NA12878.bam**
directory for temp files: /lscratch/60475214
Executing ' java -Djava.io.tmpdir=/lscratch/60475214 -Xmx2g -XX:ParallelGCThreads=2 -jar /usr/local/apps/GATK/3.8-1/GenomeAnalysisTK.jar -T CountReads -R /fdb/GATK_resource_bundle/hg38/Homo_sapiens_assembly38.fasta -I NA12878.bam --use_jdk_inflater --use_jdk_deflater '
INFO  13:50:21,043 HelpFormatter - ----------------------------------------------------------------
INFO  13:50:21,047 HelpFormatter - The Genome Analysis Toolkit (GATK) v3.8-1-0-gf15c1c3ef, Compiled...
INFO  13:50:21,048 HelpFormatter - Copyright (c) 2010-2016 The Broad Institute
INFO  13:50:21,048 HelpFormatter - For support and documentation go to https://software.broadinstit...
INFO  13:50:21,049 HelpFormatter - [Tue Jun 30 13:50:21 EDT 2020] Executing on Linux 3.10.0-862....
INFO  13:50:21,049 HelpFormatter - OpenJDK 64-Bit Server VM 1.8.0_181-b13
INFO  13:50:21,056 HelpFormatter - Program Args: -T CountReads -R /fdb/GATK_resource_bundle/hg38/Ho...
INFO  13:50:21,064 HelpFormatter - Executing as user@cn3293 on Linux 3.10.0-862.14.4.el7.x86_64 a...
INFO  13:50:21,065 HelpFormatter - Date/Time: 2020/06/30 13:50:21
INFO  13:50:21,065 HelpFormatter - ----------------------------------------------------------------
INFO  13:50:21,066 HelpFormatter - ----------------------------------------------------------------
INFO  13:50:21,072 GenomeAnalysisEngine - Deflater: JdkDeflater
INFO  13:50:21,073 GenomeAnalysisEngine - Inflater: JdkInflater
INFO  13:50:21,074 GenomeAnalysisEngine - Strictness is SILENT
INFO  13:50:23,245 GenomeAnalysisEngine - Downsampling Settings: No downsampling
INFO  13:50:23,258 SAMDataSource$SAMReaders - Initializing SAMRecords in serial
INFO  13:50:23,540 SAMDataSource$SAMReaders - Done initializing BAM readers: total time 0.28
INFO  13:50:24,955 GenomeAnalysisEngine - Preparing for traversal over 1 BAM files
INFO  13:50:24,961 GenomeAnalysisEngine - Done preparing for traversal
INFO  13:50:24,962 ProgressMeter - [INITIALIZATION COMPLETE; STARTING PROCESSING]
INFO  13:50:24,962 ProgressMeter -                 | processed |    time |    per 1M |           |   total | remaining
INFO  13:50:24,962 ProgressMeter -        Location |     reads | elapsed |     reads | completed | runtime |   runtime
INFO  13:50:24,964 ReadShardBalancer$1 - Loading BAM index data
INFO  13:50:24,966 ReadShardBalancer$1 - Done loading BAM index data
INFO  13:50:26,502 CountReads - CountReads counted 61614 reads in the traversal
INFO  13:50:26,507 ProgressMeter -            done     61614.0     1.0 s      25.0 s       99.9%     1.0 s       0.0 s
INFO  13:50:26,508 ProgressMeter - Total runtime 1.55 secs, 0.03 min, 0.00 hours
INFO  13:50:26,512 MicroScheduler - 0 reads were filtered out during the traversal out of approximately 61614 total reads (0.00%)
INFO  13:50:26,514 MicroScheduler -   -> 0 reads (0.00% of total) failing BadCigarFilter
INFO  13:50:26,515 MicroScheduler -   -> 0 reads (0.00% of total) failing MalformedReadFilter
------------------------------------------------------------------------------------------
Done. There were no warn messages.
------------------------------------------------------------------------------------------
GATK finished (exit code 0)
[user@cn3144]$ **exit**

        
```


LiftOverVCF -- converts a VCF file from one reference build to another
LiftOverVCF has been deprecated with release 3.5-0. Use picard's LiftoverVCF
 instead.


The procedure for lifting over VCF file from one genome build to different
 build in GATK is a three step process - (1) LiftoverVCF (2) sort the VCF and
 (3) FilterLiftedVCF. This process is automated by `liftOverVCF` 



```

Usage: liftOverVCF
   -vcf            <input vcf>
   -gatk           <path to gatk trunk>
   -chain          <chain file>
   -newRef         <path to new reference prefix;
                    we will need newRef.dict, .fasta, and .fasta.fai>
   -oldRef         <path to old reference prefix; we will need oldRef.fasta>
   -out            <output vcf>
   -recordOriginalLocation
                   <Should we record what the
                    original location was in the INFO field?; defaults to false>

```

This procedure is to use liftOverVCF to switch from, say, hg18 to hg19.


### Required data:


* Chain files from the
 [Broad ftp site](ftp://gsapubftp-anonymous@ftp.broadinstitute.org/Liftover_Chain_Files/)
 are mirrored in `fdb/GATK_resource_bundle/liftover_chains`. In addtion
 a hg18 to hg19 liftover file from UCSC is included as well.
* Reference sequences and .dict, .fai files for references can also be
 found in the resource bundles.


Example:



```

liftOverVCF  \
    -vcf test.hg18.vcf.gz \
    -chain /fdb/GATK_resource_bundle/liftover_chains/hg18ToHg19.over.chain \
    -newRef /fdb/GATK_resource_bundle/hg19/ucsc.hg19 \
    -oldRef /fdb/GATK_resource_bundle/hg18/Homo_sapiens_assembly18 \
    -out test.hg19.vzf

```










|  |
| --- |
| 
Quick Links
[Overview](#overview)
[Batch job on Biowulf](#serial)
[Swarm of jobs](#swarm)
[Interactive job on Biowulf](#int)
 |




Overview
### Environment module


The environment module for GATK (`module load GATK` or
 `module load GATK/version`)


* loads R as a dependency
* sets the environment variables
	+ `$GATK_TEST_DATA`


Whenever possible, please use lscratch for temporary files. It will be
 more performant and reduce strain on the shared filesystems.


### Calling GATK


GATK now includes its own wrapper script (`gatk`) and therefore
 the Biowulf-supplied wrapper for GATK 3 is not available for GATK 4. The basic 
 GATK 4 command is 



```

gatk --java-options "[java options]" <Toolname> [tool options]

```

For detailed GATK help use



```

**gatk --help**

 Usage template for all tools (uses --spark-runner LOCAL when used with a Spark tool)
    gatk AnyTool toolArgs

 Usage template for Spark tools (will NOT work on non-Spark tools)
    gatk SparkTool toolArgs  [ -- --spark-runner  sparkArgs ]

 Getting help
 gatk --list Print the list of available tools

 gatk Tool --help Print help on a particular tool

 Configuration File Specification
 --gatk-config-file PATH/TO/GATK/PROPERTIES/FILE

 gatk forwards commands to GATK and adds some sugar for submitting spark jobs

 --spark-runner  controls how spark tools are run
 valid targets are:
 LOCAL: run using the in-memory spark runner
 SPARK: run using spark-submit on an existing cluster
 --spark-master must be specified
 --spark-submit-command may be specified to control the Spark submit command
 arguments to spark-submit may optionally be specified after --
 GCS: run using Google cloud dataproc
 commands after the -- will be passed to dataproc
 --cluster  must be specified after the --
 spark properties and some common spark-submit parameters will be translated
 to dataproc equivalents

 --dry-run may be specified to output the generated command line without running it
 --java-options 'OPTION1[ OPTION2=Y ... ]' optional - pass the given string of options to the
 java JVM at runtime.
 Java options MUST be passed inside a single string with space-separated values.

```

#### Memory


The `-Xmx` option to Java specifies the maximal size of the
 heap memory. It accepts suffixes `k`, `m`, and
 `g`. This is passed to gatk in the `--java-otions` string. For
 example



```

gatk --java-options "-Xmx6g" <Toolname> [tool options]

```

Note that you should request additional memory (1-2GB) for batch/interactive jobs on
 top of the memory given to java to ensure that the whole java process stays
 below the allocated memory limit.


Note that some tools may perform better with non-default garbage collectors. For
 example here are some possible options for a large heap allocation using ParallelGCThreads:



```

gatk --java-options "-Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID -Xms48G -Xmx48G \
     -XX:ParallelGCThreads=8" ...

```

#### Redirect temporary files to /lscratch/${SLURM\_JOB\_ID}


GATK writes some temporary files during its run. By default, Java writes
 temp files to `/tmp`. However this is not advisable because
 `/tmp` is used by the operating system and is very limited in space.
 The tmp files should be redirected to `/lscratch/${SLURM_JOB_ID}`
 by adding 



```

-Djava.io.tmpdir=/lscratch/${SLURM_JOB_ID}

```

to the java options.


For more details on the node-local storage under `/lscratch` see the
 [user guide](/docs/b2-userguide.html#local)


#### Parallel garbage collection threads


If you find that your GATK runs have more active threads than you were
 expecting you may have to limit the number of prallel garbage collection
 threads. The JVM options `-XX:ParallelGCThreads` and 
 `-XX:ConcGCThreads` can be used to tune the number of threads
 dedicated to garbage collection.


### Open files


GATK may open many files simultaneously. It can happen that you hit the
 open-file limit on Biowulf, and see errors like:



```

##### ERROR MESSAGE: Couldn't read file .... (Too many open files)

```

If that happens, the limit on the number of open files has to be
 increased from its default soft limit of 4096. This can be done by adding
 the following command to your batch script: `ulimit -n 16384`.


Here are some other possibilities (thanks to Ryan Neff, NHGRI):


* Reduce the number of open files. When doing variant calling, first create gVCF
 files for each sample and then combine them together in batches of 100-200, and
 then do the final variant calls on the remaining 5-20 files (GATK's
 HaplotypeCaller and GenotypeGVCFs).



GATK Parallelization
The paradigm for parallelization for GATK 4 is completely different
 from GATK 3. Spark is used to parallelize workloads and the tools
 capable of using this form of parallelism include 'Spark' in their
 toolname. Some of these tools are still in beta and outputs may not be
 compatible with non-parallel variants. See the GATK documentation for
 more detail.




Running a single GATK batch job on Biowulf
Set up a batch script along the following lines:



```

#! /bin/bash
# gatk.bat
set -e

module load GATK/4.1.8.0 || exit 1
cd /lscratch/$SLURM_JOB_ID || exit 1

ref=/fdb/GATK_resource_bundle/hg38/Homo_sapiens_assembly38.fasta
bam=${GATK_TEST_DATA:-none}/HG00133.alt_bwamem_GRCh38DH.20150826.GBR.exome.bam
# use 8 GB memory
gatk --java-options "-Xmx8g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID" \
  CountReads \
  -R $ref \
  -I $bam \

```

Submit to the queue with [sbatch](/docs/userguide.html):



```

[user@biowulf]$ **sbatch --mem=9g --gres=lscratch:10 gatk.bat**

```

Since the simple `CountReads` tool does not support
 multithreading, only a single core (the default) is requested.




Running a swarm of GATK batch jobs on Biowulf
The following swarm commend file would run the MarkDuplicatesSpark.
 It will also sort and index the bam file. This is a Spark
 implementation of Picard MarkDuplicates that allows the tool to be run
 in parallel on multiple cores: 



```

# this file is gatk-swarm
gatk --java-options "-Xmx8g" MarkDuplicatesSpark \
  -I input1.bam \
  -O input1_markdup.bam \
  --spark-runner LOCAL \
  --spark-master local[$SLURM_CPUS_PER_TASK] \
gatk --java-options "-Xmx8g" MarkDuplicatesSpark \
  -I input2.bam \
  -O input2_markdup.bam \
  --spark-runner LOCAL \
  --spark-master local[$SLURM_CPUS_PER_TASK] \
[...]

```

These command lines specify that the process will use 8 GB of memory
 (-Xmx8g). Thus, this swarm job should be submitted with the '-g 9' flag, to
 tell swarm that each process will require 9 GB (8GB plus some buffer) of memory. 
 Please do not use more than 12 CPUs due to the parallel efficiency.



```

biowulf$ **swarm -g 9 -f gatk-swarm -t 16 --module GATK/4.1.8.0**

```



Running an interactive job on Biowulf
It may be useful for debugging purposes to run GATK jobs interactively. Such
 jobs should not be run on the Biowulf login node. Instead allocate an
 interactive node as described below, and run the interactive job there.



```

[user@biowulf]$ **sinteractive --gres=lscratch:50 --cpus-per-task=2 --mem=6g**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144]$ **cd /lscratch/$SLURM\_JOB\_ID**
[user@cn3144]$ **module load GATK/4.1.8.0**    
[user@cn3144]$ **cp ${GATK\_TEST\_DATA:-none}/NA12878.ba[im]\* .**
[user@cn3144]$ **ls -lh**
total 18M
-rw-r--r-- 1 user group 2.2M Jun 30 13:45 NA12878.bai
-rw-r--r-- 1 user group  15M Jun 30 13:45 NA12878.bam
[user@cn3144]$ **gatk --java-options "-Xmx5g -Xms5g -Djava.io.tmpdir=/lscratch/$SLURM\_JOB\_ID" \
 CountReads \
 -R /fdb/igenomes/Homo\_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa \
 -I NA12878.bam**
Using GATK jar /usr/local/apps/GATK/4.1.8.0/gatk-package-4.1.8.0-local.jar
Running:
    java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true 
    -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -Xmx5g -Xms5g 
    -Djava.io.tmpdir=/lscratch/60475214 
    -jar /usr/local/apps/GATK/4.1.8.0/gatk-package-4.1.8.0-local.jar CountReads -R 
     /fdb/igenomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa -I NA12878.bam
17:06:09.007 INFO  NativeLibraryLoader - Loading libgkl_compression.so from jar:file:/usr/local/apps/GAT
K/4.1.8.0/gatk-package-4.1.8.0-local.jar!/com/intel/gkl/native/libgkl_compression.so
17:06:09.513 INFO  CountReads - ------------------------------------------------------------
17:06:09.514 INFO  CountReads - The Genome Analysis Toolkit (GATK) v4.1.8.0
17:06:09.514 INFO  CountReads - For support and documentation go to https://software.broadinstitute.org/
gatk/
17:06:09.515 INFO  CountReads - Executing as user@cn3293 on Linux v3.10.0-862.14.4.el7.x86_64 amd64
17:06:09.515 INFO  CountReads - Java runtime: OpenJDK 64-Bit Server VM v1.8.0_181-b13
17:06:09.516 INFO  CountReads - Start Date/Time: June 30, 2020 5:06:08 PM EDT
17:06:09.517 INFO  CountReads - ------------------------------------------------------------
17:06:09.517 INFO  CountReads - ------------------------------------------------------------
17:06:09.519 INFO  CountReads - HTSJDK Version: 2.22.0
17:06:09.519 INFO  CountReads - Picard Version: 2.22.8
17:06:09.519 INFO  CountReads - HTSJDK Defaults.COMPRESSION_LEVEL : 2
17:06:09.520 INFO  CountReads - HTSJDK Defaults.USE_ASYNC_IO_READ_FOR_SAMTOOLS : false
17:06:09.520 INFO  CountReads - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_SAMTOOLS : true
17:06:09.521 INFO  CountReads - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_TRIBBLE : false
17:06:09.521 INFO  CountReads - Deflater: IntelDeflater
17:06:09.521 INFO  CountReads - Inflater: IntelInflater
17:06:09.521 INFO  CountReads - GCS max retries/reopens: 20
17:06:09.521 INFO  CountReads - Requester pays: disabled
17:06:09.521 INFO  CountReads - Initializing engine
WARNING: BAM index file /lscratch/60475214/NA12878.bai is older than BAM /lscratch/60475214/NA12878.bam
17:06:10.212 INFO  CountReads - Done initializing engine
17:06:10.212 INFO  ProgressMeter - Starting traversal
17:06:10.212 INFO  ProgressMeter -        Current Locus  Elapsed Minutes       Reads Processed     Reads
/Minute
17:06:11.192 INFO  CountReads - 0 read(s) filtered by: WellformedReadFilter

17:06:11.194 INFO  ProgressMeter - chrUn_KN707963v1_decoy:19294              0.0                 61614
      3768440.4
17:06:11.195 INFO  ProgressMeter - Traversal complete. Processed 61614 total reads in 0.0 minutes.
17:06:11.195 INFO  CountReads - CountReads counted 61614 total reads
17:06:11.195 INFO  CountReads - Shutting down engine
[June 30, 2020 5:06:11 PM EDT] org.broadinstitute.hellbender.tools.CountReads done. Elapsed time: 0.04 m
inutes.
Runtime.totalMemory()=5145362432
Tool returned:
61614

```










