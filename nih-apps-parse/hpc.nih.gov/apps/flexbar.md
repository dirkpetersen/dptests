

document.querySelector('title').textContent = ' Flexbar on Biowulf';

Flexbar on Biowulf



|  |
| --- |
| 
Quick Links
[Interactive job](#int)
[Batch job on Biowulf](#serial)
[Swarm of jobs](#swarm)
[Documentation](#doc)
 |



Description

Flexbar preprocesses high-throughput sequencing data efficiently. It 
demultiplexes barcoded runs and removes adapter sequences. Moreover, trimming 
and filtering features are provided. Flexbar increases read mapping rates and 
improves genome and transcriptome assemblies. It supports next-generation 
sequencing data in fasta and fastq format, e.g. from Illumina and the Roche 454 
platform.



### Reference


* Dodt, M., Roehr, J. T., Ahmed, R., & Dieterich, C. (2012). FLEXBAR—flexible 
 barcode and adapter processing for next-generation sequencing platforms. 
 *Biology, 1*(3), 895-905. Nicorici, D., Satalan, M., Edgren, H., 
 Kangaspeska, S., Murumagi, A.,


### Web sites


* [GitHub](https://github.com/seqan/flexbar)
* [sourceforge](http://sourceforge.net/projects/flexbar/)




Interactive job
To prepare your environment for using Flexbar, enter 
the following:




```

[user@biowulf]$ **sinteractive**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load flexbar**


```


The command **flexbar** will now be available to you. 
To learn more use the help arguments like so:




```

[user@cn3144 ~]$ **flexbar -h**

```


To copy the Flexbar test data set, use the following commands:




```

[user@cn3144 ~]$ **cp -r /usr/local/apps/flexbar/test /data/$USER/flexbar\_test\_data**
[user@cn3144 ~]$ **cd /data/$USER/flexbar\_test\_data**

```


You can now test commands with this example data set like so:




```

[user@cn3144 ~]$ **flexbar --reads test.fasta --target my\_result --adapters adapters.fasta**

               ________          __
              / ____/ /__  _  __/ /_  ____ ______
             / /_  / / _ \| |/ / __ \/ __ `/ ___/
            / __/ / /  __/>  </ /_/ / /_/ / /
           /_/   /_/\___/_/|_/_.___/\__,_/_/

Flexbar - flexible barcode and adapter removal, version 2.5


Local time:            Tue Dec 15 11:18:48 2015

Target name:           my_result
File type:             fasta
Reads file:            test.fasta
Adapter file:          adapters.fasta

threads:               1
max-uncalled:          0
min-read-length:       18

adapter-trim-end:      RIGHT
adapter-min-overlap:   3
adapter-threshold:     3
adapter-match:         1
adapter-mismatch:     -1
adapter-gap:          -6

Adapter:               Sequence:
ad1                    CGTCTT


Processing reads ...done.

Computation time:  < 1 sec


Adapter removal statistics
==========================
Adapter:            Overlap removal:    Full length:
ad1                 10                  9

Min, max, mean and median adapter overlap: 5 / 6 / 5 / 6


Output file statistics
======================
Read file:               my_result.fasta
  written reads          7
  skipped short reads    6


Filtering statistics
====================
Processed reads                   13
  skipped due to uncalled bases    0
  short prior adapter removal      2
  finally skipped short reads      6
Discarded reads overall            6
Remaining reads                    7   (53% of input)

Processed bases:   422
Remaining bases:   196   (46% of input)


Flexbar completed adapter removal.


```



Batch jobs on Biowulf

As with any sequence of commands, the last example could be written into a
script and submitted to SLURM with the sbatch command. The script may look
something like the following:




```

#!/bin/bash
# this file is called myjob.bat

module load flexbar
cd /data/$USER/flexbar_test_data
flexbar --reads test.fasta --target my_result --adapters adapters.fasta


```


Submit the job to SLURM with: 




```

[user@biowulf ~]$ **sbatch --mem=2g --walltime=00:00:10 myjob.bat**

```


Note that larger jobs may need more memory and longer walltime, and that flexbar
has a **-n** argument that controls the number of threads the program
attempts to use. If you use this argument you must ensure that your sbatch call
requests a matching number of cpus.





Swarm of jobs on Biowulf

If necessary, a swarm of Flexbar jobs could be started by writing a swarm file 
and calling swarm with the file as input. A swarm file would look something 
like so:




```

flexbar --reads 1.fasta --target result1 --adapters adapter1.fasta
flexbar --reads 2.fasta --target result2 --adapters adapter2.fasta
flexbar --reads 3.fasta --target result3 --adapters adapter3.fasta
...
flexbar --reads N.fasta --target resultN --adapters adapterN.fasta

```


Assuming this file was called **myjobs.swarm**, it would then be submitted to 
the batch system through the swarm command like so:




```

[user@biowulf ~]$ **swarm -f myjobs.swarm --module flexbar -g 2 --time 00:01:00**

```


Note that the jobs may need more memory or walltime depending on their
computational demands.
See the [swarm webpage](/apps/swarm.html) for more
information, or contact the Biowulf staff at staff@hpc.nih.gov 





Documentation
* [Flexbar manual](https://github.com/seqan/flexbar/wiki/Manual)






