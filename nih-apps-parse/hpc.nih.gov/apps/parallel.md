

document.querySelector('title').textContent = 'GNU Parallel on Biowulf';
GNU Parallel on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Important Notes](#notes)
[Examples](#examples)
 |



GNU parallel is a shell tool for executing jobs in parallel using one or more computers. A job can be a single command or a small script that has to be run for each of the lines in the input. The typical input is a list of files, a list of hosts, a list of users, a list of URLs, or a list of tables. A job can also be a command that reads from a pipe. GNU parallel can then split the input and pipe it into commands in parallel.



Documentation
* [GNU Parallel Main Site](https://www.gnu.org/software/parallel/)
* [Cheat Sheet](https://www.gnu.org/software/parallel/parallel_cheat.pdf)
* [Tutorial (2018)](https://www.gnu.org/software/parallel/parallel_tutorial.html)
* After loading the parallel module, the command

```
man parallel
```

will display a simple explanation of all options and a few examples.


* Module Name: parallel (see [the modules page](/apps/modules.html) for more information)


Important Notes
Here are some crucial options for GNU parallel:


* **-j** *N*: Run up to N jobs in parallel. Usually *N* should be set to **$SLURM\_CPUS\_PER\_TASK** when a job is submitted with either **-c** or **--cpus-per-task**.
* **-N** *max-args*: Use at most *max-args* arguments per command line.
* **--delay**  *duration*: Delay starting next job by duration.
* **--joblog** *logfile*: Logfile for executed jobs.
* **-a** or **--arg-file**: Use input-file as input source.
* **{#}**: Sequence number of the job to run.
* **--tmpdir** *dirname*: Directory for temporary files. If /lscratch is allocated, this should be set to **--tmpdir /lscratch/$SLURM\_JOB\_ID**. Otherwise, parallel will use /tmp.
* **--workdir** *dir*: Jobs will be run in the given dir, rather than the current working directory.


Examples


|  |
| --- |
| 
Quick Links
[Tar a set of directories](#ex1)
[Run a command with multiple inputs](#ex2)
[Passing multiple arguments to a parallelized command](#ex3)
[Passing large numbers of arguments](#ex4)
[Reading input from tab-separated file](#ex5)
[Saving output to unique files](#ex6)
[Combine parallel and swarm](#ex7)
 |


Here are some very useful examples. All of these examples are done after first allocating an interactive session with 8 cpus:



```
user@biowulf:~$ **sinteractive -c 8**
```

The same examples can be run in a batch job by inserting them into sbatch script:



```
#!/bin/bash
#SBATCH --cpus-per-task 8

ml parallel

...
```

### Tar a set of directories


After a project is completed, the subdirectories should be compressed into individual gzipped tar files for easy transfer and archiving. First, list the directories that are present:



```
user@node:~$ **ls -F**
001/  002/  003/  004/  005/  006/  007/  008/  README.txt  script.sh*
user@node:~$ **ls -1d 00\***
001 002 003 004 005 006 007 008
```

Tar each directory eight at a time using parallel. In this example the option **-j** tells parallel to run the given command **tar** at most eight processes simultaneously (as dictated by the **$SLURM\_CPUS\_PER\_TASK** variable set by slurm). The individual directories from **ls -1d 00\*** are passed to parallel by substituting **{}** in the command.



```
user@node:~$ **parallel -j $SLURM\_CPUS\_PER\_TASK "tar -c -z -f {}.tgz {}" ::: ls -1d 00\***
user@node:~$ **ls \*.tgz**
001.tgz  002.tgz  003.tgz  004.tgz  005.tgz  006.tgz  007.tgz  008.tgz
```

### Run a command with multiple inputs


If you have a single-threaded command that accepts an argument and outputs to a unique file, parallel can be used to accelerate the processing with a minimal of fuss.



```
user@node:~$ **parallel -j $SLURM\_CPUS\_PER\_TASK "command -i {} -o {}.out" ::: a b c d e f g h**
```

Keep in mind that because the number of simultaneous processes executed by parallel is determined by **$SLURM\_CPUS\_PER\_TASK**, the slurm option **-c** or **--cpus-per-task** will change this number. So, for the example above, if this was submitted like so:



```
sbatch -c 4 *script.sh*
```

only 4 processes would be run simultaneously.


### Passing multiple arguments to a parallelized command


The option **-N**, along with positional replacement strings, allows passing multiple arguments to a parallelized command:



```
user@node:~$ **parallel -j $SLURM\_CPUS\_PER\_TASK -N 2 "echo {1} {2}" ::: a b c d e f g h**
a b
g h
e f
c d
```

The positional replacement strings can be given in any order:



```
user@node:~$ **parallel -j $SLURM\_CPUS\_PER\_TASK -N 2 "echo {2} {1} {2}.out" ::: a b c d e f g h**
b a b.out
f e f.out
d c d.out
h g h.out
```

### Passing large numbers of arguments


If a large number of arguments needs to be passed, making the command-line unwieldy, the arguments can be written to a file. Each line is regarded as an input:



```
user@node:~$ **cat args.txt**
a
b
c
d
e
f
g
h
```

Arguments can be passed to parallel in multiple ways, either piped



```
user@node:~$ **cat args.txt | parallel -j $SLURM\_CPUS\_PER\_TASK -N 4 "echo {1} {2} {3} {4}"**
```

or by redirect



```
user@node:~$ **parallel -j $SLURM\_CPUS\_PER\_TASK -N 4 "echo {1} {2} {3} {4}" < args.txt**
```

or with the option **-a** or **--arg-file**:



```
user@node:~$ **parallel -j $SLURM\_CPUS\_PER\_TASK -N 4 -a args.txt "echo {1} {2} {3} {4}"**
```

### Reading input from tab-separated file


If a command to be parallelized has multiple arguments, it is sometimes saner to save the set of arguments as a tab-separated list in a file. For example,



```
user@node:~$ **cat args.tsv**
a       b
c       d
e       f
g       h
```


```
user@node:~$ **parallel -j $SLURM\_CPUS\_PER\_TASK -a args.tsv --colsep '\t' "echo {1}--{2}"**
a--b
c--d
g--h
e--f
```

Note that **-N** is not necessary, since **--colsep** automatically detects 2 arguments per line.


### Saving output to unique files


The **{#}** replacement string within the parallelized commands can be used to enumerate each job.



```
user@node:~$ **parallel -j $SLURM\_CPUS\_PER\_TASK -N 4 -a args.txt "echo {1} {2} {3} {4} > {#}.out"**
user@node:~$ **ls**
1.out
2.out
args.txt
user@node:~$ **cat 1.out**
a b c d
user@node:~$ **cat 2.out**
e f g h
```

### Combine parallel and swarm


[swarm](https://hpc.nih.gov/apps/swarm.html) is a tool that allows easy submission of hundreds of similiar commands in a single slurm job. In certain circumstances, those commands can be parallelized in each subjob.


Imagine a set of 80 single-threaded commands that differ in one input parameter and output is independent:



```
command -i x01.in -o 01.out
command -i x02.in -o 02.out
command -i x03.in -o 03.out
command -i x04.in -o 04.out
...
command -i x80.in -o 80.out
```

We can parallelize these to run 10 simultaneously in a swarmfile:



```
user@biowulf:~$ **cat swarmfile**
parallel -j $SLURM_CPUS_PER_TASK "command -i x{}.in -o {}.out" ::: 01 02 03 04 05 06 07 08 09 10
parallel -j $SLURM_CPUS_PER_TASK "command -i x{}.in -o {}.out" ::: 11 12 13 14 15 16 17 18 19 20
parallel -j $SLURM_CPUS_PER_TASK "command -i x{}.in -o {}.out" ::: 21 22 23 24 25 26 27 28 29 30
parallel -j $SLURM_CPUS_PER_TASK "command -i x{}.in -o {}.out" ::: 31 32 33 34 35 36 37 38 39 40
parallel -j $SLURM_CPUS_PER_TASK "command -i x{}.in -o {}.out" ::: 41 42 43 44 45 46 47 48 49 50
parallel -j $SLURM_CPUS_PER_TASK "command -i x{}.in -o {}.out" ::: 51 52 53 54 55 56 57 58 59 60
parallel -j $SLURM_CPUS_PER_TASK "command -i x{}.in -o {}.out" ::: 61 62 63 64 65 66 67 68 69 70
parallel -j $SLURM_CPUS_PER_TASK "command -i x{}.in -o {}.out" ::: 71 72 73 74 75 76 77 78 79 80
```

Then submit using swarm, allocating 10 cpus per swarm subjob:



```
user@biowulf:~$ **swarm -t 10 swarmfile**
```



