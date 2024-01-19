

document.querySelector('title').textContent = 'Swarm on Biowulf';
Swarm on Biowulf



|  |  |  |
| --- | --- | --- |
| 

|  |
| --- |
| 
Quick Links
[Video Tutorials](#videos)
[Usage](#usage)
[Details](#details)
[Input](#input)
[File Directives](#directives)
[Output](#output)
[Examples](#examples)
[STDIN/STDOUT](#stdin)
[-b, --bundle](#bundling)
[-g and -t (memory and threads)](#gandt)
[-p, --processes-per-subjob](#p)
[--time](#time)
[--dependency](#dependency)
[Fixed output path](#fixed)
[Mixed asynchronous and serial commands](#mixed)
[--module](#module)
[Setting environment variables](#environment)
[--gres (Local scratch disk space)](#lscratch)
[--sbatch](#sbatch)
[--devel, --verbose](#devel)
[Generating a swarm file](#generate)
[Monitoring a swarm](#monitor)
[Deleting/Canceling a swarm](#delete)
[Download](#download)
 |

 | 
 Swarm is a script designed to simplify submitting a group of commands to
 the Biowulf cluster. Some programs do not scale well or can't use distributed memory.
 Other programs may be 'embarrassingly parallel', in that many independent jobs need to be
 run. These programs are well suited to running 'swarms of jobs'.
 The swarm script simplifies these computational problems.
 |


Swarm reads a list of command lines (termed "commands" or "processes") from a swarm command file (termed the "swarmfile"), then automatically
submits those commands to the batch system to execute.

Command lines in the swarmfile should appear just as they would be entered on a Linux command line.

Swarm encapsulates each command line in a single temporary command script, then submits all command scripts to the Biowulf
cluster as a [Slurm job array](http://slurm.schedmd.com/job_array.html).
By default, swarm runs one command per core on a node, making optimum use of a node.
Thus, a node with 16 cores will run 16 commands **in parallel**.


For example, create a file that looks something like this (**NOTE:** lines that begin with a **#**
character are interpreted as comments and are not executed):



```
**[biowulf]$** cat file.swarm
# My first swarmfile -- this file is file.swarm
uptime
uptime
uptime
uptime
```

Then submit to the batch system:



```
**[biowulf]$** swarm --verbose 1 file.swarm
4 commands run in 4 subjobs, each command requiring 1.5 gb and 1 thread
12345
```

This will result in a single **job** (jobid 12345) of four **subjobs** (subjobids 0, 1, 2, 3), with each swarmfile line being run independently as a single subjob.
By default, each subjob is allocated a 1.5 gb of memory and 1 core (consisting of 2 cpus). 
The subjobs will be executed within the same directory from which the swarm was submitted.


The following diagram visualizes how the job array will look:



```
------------------------------------------------------------
SWARM
├── subjob 0: 1 command (1 cpu, 1.50 gb)
|   ├──  uptime
├── subjob 1: 1 command (1 cpu, 1.50 gb)
|   ├──  uptime
├── subjob 2: 1 command (1 cpu, 1.50 gb)
|   ├──  uptime
├── subjob 3: 1 command (1 cpu, 1.50 gb)
|   ├──  uptime
------------------------------------------------------------
```

All output will be written to that same directory. By default, swarm will create two output files for each independent subjob, one for 
STDOUT and one for STDERR. The format is *name*\_*jobid*\_*subjobid*.*{e,o}*:



```
**[biowulf]$** ls
file.swarm       swarm_12345_0.o  swarm_12345_1.o  swarm_12345_2.o  swarm_12345_3.o
swarm_12345_0.e  swarm_12345_1.e  swarm_12345_2.e  swarm_12345_3.e
```




Video Tutorials
[back to top](/apps/swarm.html)  




Usage
[back to top](/apps/swarm.html)  


```

Usage: swarm [swarm options] [sbatch options] swarmfile

Basic options:

  [**-g,--gb-per-process**](#gandt) [float]						       
			    gb per process (can be fractions of GB, e.g. 3.5)
  [**-t,--threads-per-process**](#gandt) [int]/"auto" 				       
			    threads per process (can be an integer or the word
			    auto).  This option is only valid for
			    multi-threaded swarms (-p 1)
  [**-p,--processes-per-subjob**](#p) [int]					       
			    processes per subjob (default = 1), this option is
			    only valid for single-threaded swarms (-t 1)
  [**-b,--bundle**](#bundling) [int]	    bundle more than one command line per subjob and
			    run sequentially (this automatically multiplies the
			    time needed per subjob)
  --noht		    don't use hyperthreading, equivalent to slurm
			    option --threads-per-core=1
  --usecsh		    use tcsh as the shell instead of bash
  --err-exit		    exit the subjob immediately on first non-zero exit
			    status
  [**-m,--module**](#module) [str]	    provide a list of environment modules to load prior
			    to execution (comma delimited)
  --no-comment		    don't ignore text following comment character #
  -c,--comment-char [char]  use something other than # as the comment character
  --maxrunning [int]	    limit the number of simultaenously running subjobs
  --merge-output	    combine STDOUT and STDERR into a single file per
			    subjob (.o)
  --logdir [dir]	    directory to which .o and .e files are to be
			    written (default is current working directory)
  --noout		    completely throw away STDOUT
  --noerr		    completely throw away STDERR
  [**--time-per-command**](#time) [str]  time per command (same as --time)
  [**--time-per-subjob**](#time) [str]   time per subjob, regardless of -b or -p

Development options:

  --no-scripts		    don't create temporary swarm scripts (with --debug
			    or --devel)
  --no-run		    don't actually run
  --debug		    don't actually run
  [**--devel**](#devel)		    combine --debug and --no-scripts, and be very
			    chatty
  [**-v,--verbose**](#devel) [int]	    can range from 0 to 6, with 6 the most verbose
  --silent		    don't give any feedback, just jobid
  -h,--help		    print this help message
  -V,--version		    print version and exit

sbatch options:

  -J,--job-name [str]	    set the name of the job
  [**--dependency**](#dependency) [str]	    set up dependency (i.e. run swarm before or after)
  [**--time**](#time) [str]		    change the walltime for each subjob (default is
			    04:00:00, or 4 hours)
  [**-L,--licenses**](/docs/userguide.html#licenses) [str]	    obtain software licenses (e.g. --licenses=matlab)
  [**--partition**](/docs/userguide.html#partitions) [str]	    change the partition (default is norm)
  [**--gres**](#lscratch) [str]		    set generic resources for swarm
  --qos [str]		    set quality of service for swarm
  --reservation [str]	    select a slurm reservation
  --exclusive		    allocate a single node per subjob, same as -t auto
  [**--sbatch**](#sbatch) [str]	    add sbatch-specific options to swarm; these options
			    will be added last, which means that swarm options
			    for allocation of cpus and memory take precedence

Environment variables:

  The following environment variables will affect how sbatch allocates
  resources:

  SBATCH_JOB_NAME           Same as --job-name
  SBATCH_TIMELIMIT          Same as --time
  SBATCH_PARTITION          Same as --partition
  SBATCH_QOS                Same as --qos
  SBATCH_RESERVATION        Same as --reservation
  SBATCH_EXCLUSIVE          Same as --exclusive

  The following environment variables are set within a swarm:

  SWARM_PROC_ID          can be 0 or 1

For more information, type "man swarm".
```




Details
[back to top](/apps/swarm.html)  



|  |  |
| --- | --- |
| swarm_fig_1 | 
A **node** consists of a hierarchy of resources.
* A **socket** is a receptacle on the motherboard for one physically packaged processor, each can contain one or more cores.
* A **core** is a complete private set of registers, execution units, and retirement queues needed to execute programs.
Nodes on the biowulf cluster can have 8, 16, or 32 cores.
* A **cpu** has the attributes of one core, but is managed and scheduled as a single logical processor by the operating system.
**Hyperthreading** is the implementation of multiple cpus on a single core.
All nodes on the biowulf cluster have hyperthreading enabled, with 2 cpus per core.
 |


Slurm allocates on the basis of **cores**. The smallest subjob runs on a single core, meaning the **smallest number of cpus that swarm can allocate is 2**.




|  |  |
| --- | --- |
| 
Swarm reads a swarmfile and creates a single **subjob** per line. By default a subjob is allocated to a single core.
Each line from a swarmfile has access to **2 cpus**.
Running swarm with the option **-t 2** is thus no different than running swarm without the -t option, as both cpus (hyperthreads)
are available to each subjob.
 | swarm_fig_2 |




|  |  |
| --- | --- |
| swarm_fig_3 | 
If commands in the swarmfile are multi-threaded, passing the -t option guarantees enough cpus will be available to the generated slurm subjobs.
For example, if the commands require either 3 or 4 threads, giving the **-t 3** or **-t 4** option allocates **2 cores per subjob**.
 |


The nodes on the biowulf cluster are configured to constrain threads within the cores the subjob is allocated. Thus, if a multi-threaded
command exceeds the cpus available, **the command will run much slower than normal!**
This may not be reflected in the overall cpu load for the node.



Memory is allocated **per subjob** by swarm, and is strictly enforced by slurm.
If a single subjob exceeds its memory allocation (by default 1.5 GB per swarmfile line), then
**the subjob will be killed by the batch system**.
See [below](#gandt) for examples on how to allocate threads and memory.




More than one swarmfile line can be run per subjob using the **-p** option. This is only valid for single-threaded
swarms (i.e. **-t 1**). Under these circumstances, all cpus are used. See [below](#p)
for more information on **-p**.






Input
[back to top](/apps/swarm.html)  

### The swarmfile


The only required argument for swarm is a swarmfile. Each line in 
the swarmfile is run as a single command. For example, the swarmfile ***file.swarm***



```
**[biowulf]$** cat file.swarm
uptime
uptime
uptime
uptime
```

when submitted like this



```
**[biowulf]$** swarm file.swarm
```

will create a swarm of 4 subjobs, with each subjob running the single command "uptime".


### Bundling


There are occasions when running a single swarmfile line per subjob is inappropriate, such as when commands
are very short (e.g. a few seconds) or when there are many thousands or millions of commands in a swarmfile. In
these circumstances, it makes more sense to ***bundle*** the swarm. For example, a swarmfile of 10,000
commands when run with a bundle value of 40 will generate 250 subjobs (10000/40 = 250):



```
**[biowulf]$** swarm --devel -b 40 file.swarm
10000 commands run in 250 subjobs, each requiring 1 gb and 1 thread, running 40 commands serially per subjob
```

**NOTE**: If a swarmfile results in more than 1000 subjobs, swarm will **automatically bundle the commands**.


**ALSO**: The time needed per subjob will be automatically multiplied by the bundle factor. If the total time
per subjob exceeds the maximum walltime of the partition, an error will be given and the swarm will not be submitted.


### Comments


By default, any text on a single line that follows a ***#*** character is assumed to be a comment,
and is ignored. For example,



```
**[biowulf]$** cat file.swarm
# Here are my commands
uptime      # this gives the current load status
pwd         # this gives the current working directory
hostname    # this gives the host name
```

However, there are some applications that require a ***#*** character in the input:



```
**[biowulf]$** cat odd.file.swarm
bogus_app -n 365#AX -w -another-flag=nonsense > output
```

The option **--no-comment** can be given to avoid removal of text following the ***#*** character.
Alternatively, another comment character can be designated using the **--comment-char** option.


### Command lists


Multiple commands can be run serially (one after the other) when they are separated by a semi-colon (;). This
is also known as a command list. For example,



```
**[biowulf]$** cat file.swarm
hostname ; date ; sleep 200 ; uptime
hostname ; date ; sleep 200 ; uptime
hostname ; date ; sleep 200 ; uptime
hostname ; date ; sleep 200 ; uptime

**[biowulf]$** swarm file.swarm
```

will create 4 subjobs, each running independently on a single cpu. Each subjob will run "hostname", followed
by "date", then "sleep 200", then "uptime", all in order.


### Complex commands


Environment variables can be set, directory locations can be changed, subshells can be spawned all within 
a single command list, and conditional statements can be given. For example, if you wanted to run some
commands in a newly created random temporary directory, you could use this:



```
**[biowulf]$** cat file.swarm
export d=/data/user/${RANDOM} ; mkdir -p $d ; if [[ -d $d ]] ; then cd $d && pwd ; else echo "FAIL" >&2 ; fi
export d=/data/user/${RANDOM} ; mkdir -p $d ; if [[ -d $d ]] ; then cd $d && pwd ; else echo "FAIL" >&2 ; fi
export d=/data/user/${RANDOM} ; mkdir -p $d ; if [[ -d $d ]] ; then cd $d && pwd ; else echo "FAIL" >&2 ; fi
export d=/data/user/${RANDOM} ; mkdir -p $d ; if [[ -d $d ]] ; then cd $d && pwd ; else echo "FAIL" >&2 ; fi
```

**NOTE**: By default, command lists are interpreted as bash commands. If a swarmfile contains tcsh- or csh-specific
commands, swarm may fail unless **--usecsh** is included.


### Line continuation markers


Application commands can be very long, with dozens of options and flags, and multiple commands separated by
semi-colons. To ease file editing, line continuation markers can be used to break up the single swarm commands
into multiple lines. For example, the swarmfile



```
cd /data/user/project; KMER="CCCTAACCCTAACCCTAA"; jellyfish count -C -m ${#KMER} -t 32 -c 7 -s 1000000000 -o /lscratch/$SLURM_JOB_ID/39sHMC_Tumor_genomic <(samtools bam2fq /data/user/bam/0A4HMC/DNA/genomic/39sHMC_genomic.md.bam ); echo ${KMER} | jellyfish query /lscratch/$SLURM_JOB_ID/39sHMC_Tumor_genomic_0 > 39sHMC_Tumor_genomic.telrpt.count
```

can be written like this:



```
cd /data/user/project; KMER="CCCTAACCCTAACCCTAA"; \
jellyfish count -C 
  -m ${#KMER} \
  -t 32 \
  -c 7 \
  -s 1000000000 \
  -o /lscratch/$SLURM_JOB_ID/39sHMC_Tumor_genomic \
  <(samtools bam2fq /data/user/bam/0A4HMC/DNA/genomic/39sHMC_genomic.md.bam ); \
echo ${KMER} | jellyfish query /lscratch/$SLURM_JOB_ID/39sHMC_Tumor_genomic_0 > 39sHMC_Tumor_genomic.telrpt.count
```

### Modules


[Environment modules](modules.html) can be loaded for an entire swarm using the **--module**
option. The 



```
swarm --module python,tophat,ucsc,samtools,vcftools -g 4 -t 8 file.swarm
```

### Swarmfile Directives


All swarm options can be incorporated into the swarmfile using swarmfile directives. Options preceded by **#SWARM** in the swarmfile (flush against the left side) will be evaluated the same as command line options.


For example, if the contents of swarmfile is as follows:



```
**[biowulf]$** cat file.swarm
#SWARM -t 4 -g 20 --time 40
command arg1
command arg2
command arg3
command arg4
```

and is submitted like so:



```
**[biowulf]$** swarm file.swarm
```

then each subjob will request 4 cpus, 20 GB of RAM and 40 minutes of walltime.


Multiple lines of swarmfile directives can be inserted, like so:



```
**[biowulf]$** cat file.swarm
#SWARM --threads-per-process 8
#SWARM --gb-per-process 8
#SWARM --sbatch '--mail-type=FAIL --export=var=100,nctype=12 --chdir=/data/user/test'
#SWARM --logdir /data/user/swarmlogs
command
command
command
command
```

The precedence for options is handled in the same way as sbatch, but with options provided with the **--sbatch** option last:



```
    command line > environment variables > swarmfile directives > --sbatch options 
```

Thus, if the swarmfile has:



```
**[biowulf]$** cat file.swarm
#SWARM -t 4 -g 20 --time 40 --partition norm
command arg1
command arg2
command arg3
command arg4
```

and is submitted like so:



```
**[biowulf]$** SBATCH_PARTITION=quick swarm -g 10 --time=10 file.swarm
```

then each subjob will request 4 cpus, 10 GB of RAM and 10 minutes of walltime. The amount of memory and walltime requested with command line options and the partition chosen with the **SBATCH\_PARTITION** environment variable supersedes the amount requested with swarmfile directives.


**NOTE:** All lines with correctly formatted **#SWARM** directives will be removed even if **--no-comment** or a non-default **--comment-char** is given.





Output
[back to top](/apps/swarm.html)  

### Default output files


STDOUT and STDERR output from subjobs executed under swarm will be
directed to a file named **swarm\_*jobid\_subjobid*.o** and **swarm\_*jobid\_subjobid*.e**, respectively. 


Please pay attention to the memory requirements of your swarm jobs!

When a swarm job runs out of memory, the node stalls and the job is eventually killed or
dies.

At the bottom of the .e file, you may see a warning like this:



```
slurmstepd: Exceeded job memory limit at some point. Job may have been partially swapped out to disk.
```

 If a job dies before it is finished, this output may not be available. Contact
[staff@hpc.nih.gov](mailto:staff@hpc.nih.gov) when you have a question about why
a swarm stopped prematurely.


### Renaming output files


The sbatch option **--job-name** can be used to rename the default output files.



```
**[biowulf]$** swarm -f file.swarm --job-name programAOK
...
**[biowulf]$** ls
programAOK_21381_0.e  programAOK_21381_2.e  programAOK_21381_4.e  programAOK_21381_6.e
programAOK_21381_0.o  programAOK_21381_2.o  programAOK_21381_4.o  programAOK_21381_6.o
programAOK_21381_1.e  programAOK_21381_3.e  programAOK_21381_5.e  programAOK_21381_7.e
programAOK_21381_1.o  programAOK_21381_3.o  programAOK_21381_5.o  programAOK_21381_7.o
```

### Combining STDOUT and STDERR into a single file per subjob


Including the **--merge-output** option will cause the STDERR output to be combined into the file used
for STDOUT. For swarm, that means the content of the .e files are written to the .o file. Keep in mind that
interweaving of content will occur.



```
**[biowulf]$** swarm --merge-output file.swarm
...
**[biowulf]$** ls
swarm_50158339_0.o   swarm_50158339_1.o  swarm_50158339_4.o  swarm_50158339_7.o
swarm_50158339_10.o  swarm_50158339_2.o  swarm_50158339_5.o  swarm_50158339_8.o
swarm_50158339_11.o  swarm_50158339_3.o  swarm_50158339_6.o  swarm_50158339_9.o
```

### Writing output files to a separate directory


By default, the STDOUT and STDERR files are written to the same directory from which the swarm
was submitted. To redirect the files to a different directory, use **--logdir**:



```
swarm --logdir /path/to/another/directory file.swarm
```

### Redirecting output


Input/output redirects (and everything in the swarmfile) should be bash compatible. For example,



```
**[biowulf]$** cat bash_file.swarm
program1 -o -f -a -n 1 > output1.txt 2>&1
program1 -o -f -a -n 2 > output2.txt 2>&1
**[biowulf]$** swarm bash_file.swarm
```

csh-style redirects like '**program >&; output**' will not work correctly unless
the **--usecsh** option is included. For example,



```
**[biowulf]$** cat csh_file.swarm
program1 -o -f -a -n 1 >& output1.txt
program1 -o -f -a -n 2 >& output2.txt
**[biowulf]$** swarm **--usecsh** csh_file.swarm
```

Be aware of programs that write directly to a file using a fixed filename.
A file will be overwritten and garbled if multiple processes are writing to the same file.
If you run multiple instances of such programs then for each instance you will
need to a) change the name of the file in each command **or** b) alter the path to the file. See
the **EXAMPLES** section for some ideas.





Examples
[back to top](/apps/swarm.html)  



|  |  |  |
| --- | --- | --- |
| 

|  |
| --- |
| 
Quick Links
[STDIN/STDOUT](#stdin)
[-b, --bundle](#bundling)
[-g and -t](#gandt)
[-p, --processes-per-subjob](#p)
[--time](#time)
[--dependency](#dependency)
[Fixed output path](#fixed)
[Mixed asynchronous and serial commands](#mixed)
[--module](#module)
[Setting environment variables](#environment)
[--sbatch](#sbatch)
[--devel, --verbose](#devel)
 |

 | 
 To see how swarm works, first create a file containing a few simple
 commands, then use swarm to submit them to the batch queue:
 |



```
**[biowulf]$** cat > file.swarm
date
hostname
ls -l
^D
**[biowulf]$** swarm file.swarm
```

Use **sjobs** to monitor the status of your request; an
"R" in the "St"atus column indicates your job is running.
This particular example will probably run to completion
before you can give the qstat command. To see the output from the commands, see
the files named **swarm\_*#\_#*.o**.




---


[back to top](/apps/swarm.html)





**A program that reads to STDIN and writes to STDOUT**


For each invocation of the program the names for the input and output files
vary:



```
**[biowulf]$** cat > runbix
./bix < testin1 > testout1
./bix < testin2 > testout2
./bix < testin3 > testout3
./bix < testin4 > testout4
^D
```



---


[back to top](/apps/swarm.html)



**Bundling large numbers of commands**


If you have over 1000 commands, especially if each one runs for a short
time, you should 'bundle' your jobs with the **-b** flag. For example, if the 
swarmfile contains 2560 commands, the following swarm command will group them into
bundles of 40 commands each, producing 64 command bundles. Swarm will then submit the 
64 command bundles, rather than the 2560 commands individually, as a single swarm job.
This would result in a swarm of 64 (2560/40) subjobs.



```
**[biowulf]$** swarm -b 40 file.swarm
```

Note that commands in a bundle will run sequentially on the assigned node.



---


[back to top](/apps/swarm.html)





**Allocating memory and threads with -g and -t options**


If the subjobs require significant amounts of memory (> 1.5 GB) or threads (> 1 per core), a swarm can
run fewer subjobs per node than the number of cores available
on a node. For example, if the commands in a swarmfile need up to 40 GB of
memory each using 8 threads, running swarm with --devel shows what might happen:



```
**[biowulf]$** swarm -g 40 -t 8 --devel file.swarm
14 commands run in 14 subjobs, each requiring 40 gb and 8 threads
```

If a command requires to use as many cpus on a node as possible, then the option **-t auto**
should be added. This causes each subjob in the swarmfile to allocate an entire node exclusively to the subjob, allowing
the subjob to use all available cpus on the node.


The default partition **norm** has nodes with a maximum of 248GB memory. If **-g** exceeds 373GB, swarm will give a warning message:



```
**[biowulf]$** swarm -g 400 file.swarm
ERROR: -g 400 requires --partition largemem
```

To allocate more than 373GB of memory per command, include **--partition largmem**:



```
**[biowulf]$** swarm -g 500 --partition largemem file.swarm
```

For more information about partitions, please see [https://hpc.nih.gov/docs/userguide.html#partitions](/docs/userguide.html#partitions)




---


[back to top](/apps/swarm.html)



**Using -p option to "pack" commands**


By default, swarm allocates a single command line per subjob. If the command is single-threaded, then swarm wastes half the
cpus allocated, because the slurm batch system allocates no less than a single core (or two cpus) per subjob. This effect
can be seen using the **jobload** command for a 4-command swarm:



```
**[biowulf]$** swarm file.swarm
219433
**[biowulf]$**$ jobload -u user
         JOBID         TIME   NODES   CPUS  THREADS   LOAD             MEMORY
                                     Alloc  Running                Used/Alloc
      219433_3         0:37  cn0070      2        1    50%      1.0 GB/1.5 GB
      219433_2         0:37  cn0070      2        1    50%      1.0 GB/1.5 GB
      219433_1         0:37  cn0069      2        1    50%      1.0 GB/1.5 GB
      219433_0         0:37  cn0069      2        1    50%      1.0 GB/1.5 GB

USER SUMMARY
     Jobs: 2
    Nodes: 2
     CPUs: 4
 Load Avg: 50%
```

In order to use all the cpus allocated to a single-threaded swarm, the option **-p** will set the number
of commands run per subjob. Including **-p 2**, half as many subjobs are created, each using twice as many cpus and twice as much memory:



```
**[biowulf]$** swarm -p 2 file.swarm
219434
**[biowulf]$**$ jobload -u user
         JOBID         TIME   NODES   CPUS  THREADS   LOAD             MEMORY
                                     Alloc  Running                Used/Alloc
      219434_1         0:24  cn0069      2        2   100%      2.0 GB/3.0 GB
      219434_0         0:24  cn0069      2        2   100%      2.0 GB/3.0 GB

USER SUMMARY
     Jobs: 2
    Nodes: 2
     CPUs: 4
 Load Avg: 100%
```

In this case, we are "packing" 2 commands per subjob.


**NOTE:** The cpus on the biowulf cluster are *hypercores*, and some programs run more inefficiently
when packed onto hypercores. Please test your application to see if it actually benefits from running two commands per core rather than one.


Keep in mind:


* **-p** is only available to single-threaded swarms (i.e. **-t 1**).
* The default file output format is different using **-p**. The file names end with an extra suffix indicating the cpu from the subjob:



```

**[biowulf]$**$ swarm -p 2 ../file.swarm
14 commands run in 7 subjobs, each command requiring 1.5 gb and 1 thread, packing 2 processes per subjob
221574
**[biowulf]$**$ ls
swarm_221574_0_0.e  swarm_221574_1_1.e  swarm_221574_3_0.e  swarm_221574_4_1.e  swarm_221574_6_0.e
swarm_221574_0_0.o  swarm_221574_1_1.o  swarm_221574_3_0.o  swarm_221574_4_1.o  swarm_221574_6_0.o
swarm_221574_0_1.e  swarm_221574_2_0.e  swarm_221574_3_1.e  swarm_221574_5_0.e  swarm_221574_6_1.e
swarm_221574_0_1.o  swarm_221574_2_0.o  swarm_221574_3_1.o  swarm_221574_5_0.o  swarm_221574_6_1.o
swarm_221574_1_0.e  swarm_221574_2_1.e  swarm_221574_4_0.e  swarm_221574_5_1.e
swarm_221574_1_0.o  swarm_221574_2_1.o  swarm_221574_4_0.o  swarm_221574_5_1.o
```

In the case where each swarm subjob must create or use a unique directory or file, an environment variable **SWARM\_PROC\_IC** is
available to discriminate the 0 and 1 processes running with -p 2.


For example, in order to create a unique directory in allocated /lscratch for each subjob, this bash code example can be used:



```

export TAG=${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}_${SWARM_PROC_ID} && mkdir /lscratch/${SLURM_JOB_ID}/${TAG} && touch /lscratch/${SLURM_JOB_ID}/${TAG}/foo.{0..4} && tar czf /data/user/${TAG}.tgz /lscratch/${SLURM_JOB_ID}/${TAG}/foo.*
export TAG=${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}_${SWARM_PROC_ID} && mkdir /lscratch/${SLURM_JOB_ID}/${TAG} && touch /lscratch/${SLURM_JOB_ID}/${TAG}/foo.{0..4} && tar czf /data/user/${TAG}.tgz /lscratch/${SLURM_JOB_ID}/${TAG}/foo.*
export TAG=${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}_${SWARM_PROC_ID} && mkdir /lscratch/${SLURM_JOB_ID}/${TAG} && touch /lscratch/${SLURM_JOB_ID}/${TAG}/foo.{0..4} && tar czf /data/user/${TAG}.tgz /lscratch/${SLURM_JOB_ID}/${TAG}/foo.*
export TAG=${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}_${SWARM_PROC_ID} && mkdir /lscratch/${SLURM_JOB_ID}/${TAG} && touch /lscratch/${SLURM_JOB_ID}/${TAG}/foo.{0..4} && tar czf /data/user/${TAG}.tgz /lscratch/${SLURM_JOB_ID}/${TAG}/foo.*
export TAG=${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}_${SWARM_PROC_ID} && mkdir /lscratch/${SLURM_JOB_ID}/${TAG} && touch /lscratch/${SLURM_JOB_ID}/${TAG}/foo.{0..4} && tar czf /data/user/${TAG}.tgz /lscratch/${SLURM_JOB_ID}/${TAG}/foo.*
export TAG=${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}_${SWARM_PROC_ID} && mkdir /lscratch/${SLURM_JOB_ID}/${TAG} && touch /lscratch/${SLURM_JOB_ID}/${TAG}/foo.{0..4} && tar czf /data/user/${TAG}.tgz /lscratch/${SLURM_JOB_ID}/${TAG}/foo.*
export TAG=${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}_${SWARM_PROC_ID} && mkdir /lscratch/${SLURM_JOB_ID}/${TAG} && touch /lscratch/${SLURM_JOB_ID}/${TAG}/foo.{0..4} && tar czf /data/user/${TAG}.tgz /lscratch/${SLURM_JOB_ID}/${TAG}/foo.*

```

In this case, while the files created within each distinct **/lscratch/${SLURM\_JOB\_ID}/${SLURM\_ARRAY\_JOB\_ID}\_${SLURM\_ARRAY\_TASK\_ID}\_${SWARM\_PROC\_ID}** directory are identical to all
the other swarm subjobs, the final tarball is unique:



```

**[biowulf]$**$ ls /data/user/output
221574_0_0.tgz  221574_0_1.tgz  221574_1_0.tgz  221574_1_1.tgz  221574_2_0.tgz  221574_2_1.tgz
221574_3_0.tgz  221574_3_1.tgz  221574_4_0.tgz  221574_4_1.tgz  221574_5_0.tgz  221574_5_1.tgz 
221574_6_0.tgz  221574_6_1.tgz

```



---


[back to top](/apps/swarm.html)



**Setting walltime with --time**


By default all jobs and subjobs have a walltime of 2 hours. If a swarm subjob exceeds its walltime, **it will be killed!**.
On the other hand, if your swarm subjobs have a very short walltime, then their priority on the queue may be elevated. Therefore,
it is best practice to set a walltime using the **--time** option that reflects the estimated execution time of the subjobs.
For example, if the command lines in a swarm are expected to require no more than half an hour to complete, the swarm command should be:



```
**[biowulf]$** swarm --time 00:30:00 file.swarm
```

Because a subjob is expected to be running a single command from the swarmfile, the value of **--time** can be considered
the amount of time to run a single command. When a swarm is bundled, the value for **--time** is then
multiplied by the bundle factor. For example, if
a swarm that normally creates 64 commands is bundled to run 4 commands serially, the value of **--time** is
multiplied by 4:



```
**[biowulf]$** swarm **--time 00:30:00 -b 4** --devel file.swarm
64 commands run in 16 subjobs, each command requiring 1.5 gb and 1 thread, running 4 processes serially per subjob
sbatch --array=0-15 --job-name=swarm **--time=2:00:00** --cpus-per-task=2 --partition=norm --mem=1536
```

If a swarm has more than 1000 commands and is autobundled, there is a chance that the time requested will exceed
the maximum allowed. In that case, an error will be thrown:



```

ERROR: Total time for bundled commands is greater than partition walltime limit.
Try lowering the time per command (--time=04:00:00), lowering the bundle factor
(if not autobundled), picking another partition, or splitting up the swarmfile.
```

See the  [Biowulf User Guide for a discussion of walltime limits](/docs/userguide.html#wall).


There are two additional options for setting the time of a swarm. **--time-per-command** is identical to **--time**, and 
merely serves as a more obvious explanation of time allocation.


**--time-per-subjob** overrides the time adjustments applied when [bundling](#bundling) or [packing](#p) commands.
This option can be used when a single command takes less than 1 minute to complete and there are a high number of commands bundled per
subjob:



```
**[biowulf]$** swarm **--time-per-subjob 00:30:00 -b 4** --devel file.swarm
64 commands run in 16 subjobs, each command requiring 1.5 gb and 1 thread, running 4 processes serially per subjob
sbatch --array=0-15 --job-name=swarm **--time=30:00** --cpus-per-task=2 --partition=norm --mem=1536
```



---


[back to top](/apps/swarm.html)



**Handling job dependencies**



If a swarm is run as a single step in a pipeline, job dependencies can be handled with the **--dependency** options.
For example, a first script (first.sh) is to be run to generate some initial data files. Once this job is finished, a swarm of
commands (swarmfile.txt) is run to take the output of the first script and process it. Then, a last script (last.sh) is run
to consolidate the output of the swarm and further process it into its final form.




Below, the swarm is run with a dependency on the first script. Then the last script is run with a dependency on the swarm.
The swarm will sit in a pending state until the first job (10001) is completed, and the last job will sit in a pending state until
the entire swarm (10002) is completed.




```
**[biowulf]$** sbatch first.sh
10001
**[biowulf]$** swarm --dependency afterany:10001 file.swarm
10002
**[biowulf]$** sbatch --dependency=afterany:10002 last.sh
10003
```


The jobid of a job can be captured from the sbatch command and passed to subsequent submissions in a script (master.sh).
For example, here is a bash script which automates the above procedure, passing the variable $id to the first script. In this way,
the master script can be reused for different inputs:


```
**[biowulf]$** cat master.sh
#!/bin/bash
jobid1=$(sbatch first.sh)
echo $jobid1
jobid2=$(swarm --dependency afterany:$jobid1 file.swarm)
echo $jobid2
jobid3=$(sbatch --dependency=afterany:$jobid2 last.sh)
echo $jobid3
```

Now, master.sh can be submitted with a single argument



```
**[biowulf]$** bash master.sh mydata123
10001
10002
10003
**[biowulf]$**
```

You can check on the job status using squeue:

```
**[biowulf]$** squeue -u user
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
       10002_[0-3]      norm    swarm     user PD       0:00      1 (Dependency)
             10003      norm  last.sh     user PD       0:00      1 (Dependency)
             10001      norm first.sh     uwer  R       0:33      1 cn0121
```


The dependency key 'afterany' means run only after the entire job finishes, regardless of its exit status. Swarm passes the exit
status of the last command executed back to Slurm, and Slurm consolidates all the exit statuses of the subjobs in the job array into
a single exit status.



The final statuses for the jobs can be seen with sacct. The individual subjobs from swarm are designated
by **jobid\_subjobid**:



```
**[biowulf]$** sacct
       JobID    JobName  Partition    Account  AllocCPUS      State ExitCode
------------ ---------- ---------- ---------- ---------- ---------- --------
10001          first.sh       norm       user          2  COMPLETED      0:0
10001.batch       batch                  user          1  COMPLETED      0:0
10002_3           swarm       norm       user          2     FAILED      2:0
10002_3.bat+      batch                  user          1     FAILED      2:0
10003           last.sh       norm       user          2  COMPLETED      0:0
10003.batch       batch                  user          1  COMPLETED      0:0
10002_0           swarm       norm       user          2  COMPLETED      0:0
10002_0.bat+      batch                  user          1  COMPLETED      0:0
10002_1           swarm       norm       user          2  COMPLETED      0:0
10002_1.bat+      batch                  user          1  COMPLETED      0:0
10002_2           swarm       norm       user          2  COMPLETED      0:0
10002_2.bat+      batch                  user          1  COMPLETED      0:0
```

If any of the subjobs in the swarm failed, the job is marked as **FAILED**. It almost all cases, it is better
to rely on **afterany** rather than **afterok**, since the latter may cause the dependent job to 
remain queued forever:



```
**[biowulf]$** sjobs
                                                       ................Requested............................
User       JobId    JobName    Part  St      Runtime   Nodes  CPUs     Mem        Dependency     Features             Nodelist
user       10003    last.sh    norm   PD          0:00   1       1   2.0GB/cpu   afterok:10002_*   (null)               (DependencyNeverSatisfied)
```

See [the Biowulf User Guide](/docs/userguide.html#depend), or [SchedMD for a discussion on how Slurm handles exit codes](http://slurm.schedmd.com/job_exit_code.html).


NOTE: Setting **-p** causes multiple commands to run per subjob. Because of this, the exit status of the
subjob can come from any of the multiple processes in the subjob. 




---


[back to top](/apps/swarm.html)





**A program that writes to a fixed filepath**


If a program writes to a fixed filename, then you may need to run the
program in different directories. First create the necessary directories (for
instance run1, run2), and in the swarmfile cd to the unique output
directory before running the program: (cd using either an absolute path
beginning with "/" or a relative path from your home directory). Lines with
leading "#" are considered comments and ignored.



```
**[biowulf]$** cat > file.swarm
# Run ped program using different directory
# for each run
cd pedsystem/run1; ../ped
cd pedsystem/run2; ../ped
cd pedsystem/run3; ../ped
cd pedsystem/run4; ../ped
...

**[biowulf]$** swarm file.swarm
```



---


[back to top](/apps/swarm.html)





**Running mixed asynchronous and serial commands in a swarm**


There are occasions when a single swarm command can contain a mixture of asynchronous and serial commands. For 
example, collating the results of several commands into a single output and then running another command on the pooled
results. If run interactively, it would look like this:



```

**[biowulf]$** cmdA < inp.1 > out.1
**[biowulf]$** cmdA < inp.2 > out.2
**[biowulf]$** cmdA < inp.3 > out.3
**[biowulf]$** cmdA < inp.4 > out.4
**[biowulf]$** cmdB -i out.1 -i out.2 -i out.3 -i out.4 > final_result

```

It would be more efficient if the four **cmdA** commands could run asynchronously (in parallel), and then
the last **cmdB** command would wait until they were all done and then run, all on the same node and in the same
swarm command. This can be achieved using process substitution with this one-liner in a swarmfile:



```

( cmdA < inp.1 > out.1 & cmdA < inp.2 > out.2 & \
  cmdA < inp.3 > out.3 & cmdA < inp.4 > out.4 & wait ) ; \
  cmdB -i out.1 -i out.2 -i out.3 -i out.4 > final_result

```

Here, the **cmdA** commands are all run asynchronously in four background processes, and the **wait** command
is given to prevent **cmdB** from running until all the background processes are finished. Note that line
continuation markers were used for easier editing.




---


[back to top](/apps/swarm.html)



**Using --module option**


It is sometimes difficult to set the environment properly before running commands. The
easiest way to do this on Biowulf is with [environment modules](modules.html). Running commands via swarm complicates the issue, because the modules
must be loaded prior to every line in the swarmfile. Instead, you can use the **--module**
option to load a list of modules:



```
**[biowulf]$** swarm --module ucsc,matlab,python/2.7 file.swarm
```

Here, the environment is set to use the UCSC executables, Matlab, and an older, non-default
version of Python.




---


[back to top](/apps/swarm.html)



**Using --gres option**


[Local scratch disk space is NOT automatically available under Slurm](http://hpc.nih.gov/docs/userguide.html#local). Instead, local scratch disk space
is allocated using **--gres**. Here is an example of how to allocate 200GB of local scratch disk space for each swarm command:



```
**[biowulf$** swarm --gres=lscratch:200 file.swarm
```

Including **--gres=lscratch:*N***, where ***N*** is the number of GB required, will create a subdirectory on the node
corresponding to the jobid, e.g.:



```
**/lscratch/987654/**
```

This local scratch directory can be accessed dynamically using the **$SLURM\_JOB\_ID** environment variable:



```
**/lscratch/$SLURM\_JOB\_ID/**
```

Local scratch space is allocated **per subjob**. By default, that means each command or command list (single line in
swarmfile) is allocated its own independent local scratch space. **HOWEVER**, there are two situations where some
thought must be given to local scratch space:


* **bundled swarms** - [Bundled swarms](#bundling) serialize multiple commands into a single subjob. Since local scratch space is
allocated per subjob, this means that each command in the job inherits the same local scratch space. This means that each
command should be written to deal with any "leftover" files from the previous commands. A simple solution might be to
clean out the local scratch space at the end of each command. For example:  


```
cd /lscratch/$SLURM_JOB_ID ; command1 arg1 arg2 ; rm -rf /lscratch/$SLURM_JOB_ID/*
```
* **-p 2** - If the **[-p 2](#p)** option is given to swarm, then allocated local scratch space is shared
between 2 commands in a single job. In this case, make sure to allocate **twice** as much local scratch space as
normal.




---


[back to top](/apps/swarm.html)



**Setting environment variables**


If an entire swarm requires one or more environment variables to be set, the sbatch option **--export**
can be used to set the variables prior to running. In this example, we need to set the BOWTIE\_INDEXES environment variable
to the correct path for all subjobs in the swarm:



```
**[biowulf]$** swarm --sbatch "--export=BOWTIE_INDEXES=/fdb/igenomes/Mus_musculus/UCSC/mm9/Sequence/BowtieIndex/" file.swarm
```

**NOTE:** Environment variables set with the **--sbatch "--export="** option are defined
**PRIOR** to the job being submitted. This prevents setting environment variables using Slurm-generated environment
variables, such as $SLURM\_JOB\_ID or $SLURM\_MEM\_PER\_NODE.


However, if each command line in the swarm requires a unique set of environment variables, this must be done in the swarmfile. For example, setting TMPDIR to a unique subdirectory of /lscratch/$SLURM\_JOB\_ID:



```
**[biowulf]$** cat file.swarm 
export TMPDIR=/lscratch/$SLURM_JOB_ID/xyz1; mkdir $TMPDIR; cmdxyz -x 1 -y 1 -z 1
export TMPDIR=/lscratch/$SLURM_JOB_ID/xyz2; mkdir $TMPDIR; cmdxyz -x 2 -y 2 -z 2
export TMPDIR=/lscratch/$SLURM_JOB_ID/xyz3; mkdir $TMPDIR; cmdxyz -x 3 -y 3 -z 3
export TMPDIR=/lscratch/$SLURM_JOB_ID/xyz4; mkdir $TMPDIR; cmdxyz -x 4 -y 4 -z 4
```

Further, if individual commands within each command line require unique environment variables, this can be done by prefacing the command itself with the variable set:



```
**[biowulf]$** cat file.swarm
MYENV=1 command ; MYENV=2 command ; MYENV=3 command
MYENV=4 command ; MYENV=5 command ; MYENV=6 command
```



---


[back to top](/apps/swarm.html)



**Using sbatch flags**


Swarm creates a single [job array via the sbatch command](http://hpc.nih.gov/docs/userguide.html#submit)
; all valid sbatch commandline options are also valid for swarm.
However, they must be passed with a single **--sbatch** option, surrounded
by quotation marks. In this example some extra sbatch options are added.



```
**[biowulf]$** swarm --sbatch "--mail-type=FAIL --export=var=100,nctype=12 --chdir=/data/user/test" file.swarm
```

In this case,
 * **--mail-type=FAIL**: causes a single email per swarm to be sent if one subjob fails for some reason
* **--mail-type=END,ARRAY\_TASKS**: causes each subjob to issue an email when the subjob ends for any reason
Unless the ARRAY\_TASKS option is specified, mail notifications on job BEGIN, END and FAIL apply to a job array as a whole rather than generating individual email messages for each task in the job array. See [userguide.html#email](/docs/userguide.html#email) for more information about email notifications.
* **--export=var=100,nctype=12**: sets two environment variables ***var*** and ***nctype*** to 100 and 12 prior to running
* **--chdir=/data/user/test**: relocate to that directory prior to running any commands





**NOTE:** Sbatch options passed through swarm using the **--sbatch** option are listed last in the sbatch command, and thus will **override** swarm
options. When [bundling](#bundling) or [packing](#p) commands, **DO NOT** use **--time, --cpus-per-task, --mem, or --mem-per-cpu** sbatch options, as they will
inevitably conflict with those values set by swarm per command.




---


[back to top](/apps/swarm.html)



**Using --devel and --verbose options**


Before submitting a large complex swarm to the batch system, it is better
to see what would happen before it's too late. In this case, the **--devel** option
will display a good deal of information. This example shows a **huge** number of
commands autobundled to run 346 command lines serially per core.



```
**[biowulf]$** swarm --devel file.swarm
345029 commands run in 998 subjobs, each requiring 1 gb and 1 thread, running 346 commands serially per subjob
```

**--verbose** accepts a number between 0 (the same as **--silent**) and 6. Increasing the verbosity level with **--verbose** and including **--devel** will give a visual representation
of the swarm, along with lots of information about the swarm:



```
**[biowulf]$** swarm --devel --verbose 5 -g 5 -p 4 -b 4 file.swarm
basedir = /spin1/swarm/user
script dir = /spin1/swarm/user/cF2_0V7N
------------------------------------------------------------
SWARM
├── subjob 0: 16 commands (4 cpus, 20.00 gb)
|   ├── cmd 0 ; cmd 1 ; cmd 2 ; cmd 3 ;
|   ├── cmd 4 ; cmd 5 ; cmd 6 ; cmd 7 ;
|   ├── cmd 8 ; cmd 9 ; cmd 10 ; cmd 11 ;
|   ├── cmd 12 ; cmd 13 ; cmd 14 ; cmd 15 ;
├── subjob 1: 16 commands (4 cpus, 20.00 gb)
|   ├── cmd 16 ; cmd 17 ; cmd 18 ; cmd 19 ;
|   ├── cmd 20 ; cmd 21 ; cmd 22 ; cmd 23 ;
|   ├── cmd 24 ; cmd 25 ; cmd 26 ; cmd 27 ;
|   ├── cmd 28 ; cmd 29 ; cmd 30 ; cmd 31 ;
------------------------------------------------------------
2 subjobs, 32 commands, 8 output files
32 commands run in 2 subjobs, each command requiring 5 gb and 1 thread, packing 4 processes per subjob, running 4 processes serially per subjob
sbatch --array=0-1 --job-name=swarm --output=/dev/null --error=/dev/null --cpus-per-task=4 --mem=20480 /spin1/swarm/user/cF2_0V7N.batch
```

This shows a swarm of 32 commands (show as "cmd 0" ==> "cmd 31") within 2 subjobs. Each command requires 5 gb of memory,
and the commands are bundled to run 4 commands sequentially on the cpus allocated.




---


[back to top](/apps/swarm.html)



Generating a swarm script
[back to top](/apps/swarm.html)  
