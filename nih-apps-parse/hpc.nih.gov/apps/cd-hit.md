

document.querySelector('title').textContent = 'cd-hit on Biowulf';
cd-hit on Biowulf


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



CD-HIT is a very widely used program for clustering and comparing protein or nucleotide sequences. CD-HIT is very fast and can handle extremely large databases. CD-HIT helps to significantly reduce the computational and manual efforts in many sequence analysis tasks and aids in understanding the data structure and correct the bias within a dataset. 



### References:


* CD-HIT was originally developed
by [Dr. Weizhong Li](mailto:liwz@sdsc.edu) at
[Dr. Adam Godzik's Lab](http://bioinformatics.burnham.org/) at 
[the Burnham Institute (now Sanford-Burnham Medical Research Institute)](http://www.sanfordburnham.org/)


Documentation
* [cd-hit website](http://weizhongli-lab.org/cd-hit/)
* [Documentation at CD-HIT wiki](https://github.com/weizhongli/cdhit/wiki)


Important Notes
* Module Name: cd-hit (see [the modules page](/apps/modules.html) for more information)
* Multithreaded



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive --cpus-per-task=4**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load cd-hit**

[user@cn3144 ~]$ **cd-hit -i /fdb/fastadb/drosoph.aa.fas -o drosoph100 -c 1.00 -n 5 -M 16000 -d 0 -T $SLURM\_CPUS\_PER\_TASK**
================================================================
Program: CD-HIT, V4.7 (+OpenMP), May 10 2018, 13:51:01
Command: cd-hit -i /fdb/fastadb/drosoph.aa.fas -o drosoph100
         -c 1.00 -n 5 -M 16000 -d 0 -T 4

Started: Thu May 10 13:55:50 2018
================================================================
                            Output
----------------------------------------------------------------
Warning: total number of CPUs in the system is 2
Actual number of CPUs to be used: 2

total seq: 14329
longest and shortest : 8805 and 11
Total letters: 7178839
Sequences have been sorted

Approximated minimal memory consumption:
Sequence        : 9M
Buffer          : 2 X 12M = 25M
Table           : 2 X 65M = 131M
Miscellaneous   : 0M
Total           : 165M

Table limit with the given memory limit:
Max number of representatives: 4000000
Max number of word counting entries: 1979299596

# comparing sequences from          0  to       3582
...---------- new table with     3409 representatives
# comparing sequences from       3582  to       6268
99.9%---------- new table with     2559 representatives
# comparing sequences from       6268  to       8283
----------    265 remaining sequences to the next cycle
---------- new table with     1680 representatives
# comparing sequences from       8018  to       9595
----------    317 remaining sequences to the next cycle
---------- new table with     1215 representatives
# comparing sequences from       9278  to      10540
..........    10000  finished       9566  clusters
----------    446 remaining sequences to the next cycle
---------- new table with      794 representatives
# comparing sequences from      10094  to      11152
----------    269 remaining sequences to the next cycle
---------- new table with      769 representatives
# comparing sequences from      10883  to      11744
----------    207 remaining sequences to the next cycle
---------- new table with      625 representatives
# comparing sequences from      11537  to      12235
----------    169 remaining sequences to the next cycle
---------- new table with      516 representatives
# comparing sequences from      12066  to      12631
----------    121 remaining sequences to the next cycle
---------- new table with      430 representatives
# comparing sequences from      12510  to      12964
----------    103 remaining sequences to the next cycle
---------- new table with      346 representatives
# comparing sequences from      12861  to      13228
----------     76 remaining sequences to the next cycle
---------- new table with      284 representatives
# comparing sequences from      13152  to      13446
----------     60 remaining sequences to the next cycle
---------- new table with      227 representatives
# comparing sequences from      13386  to      14329
.....................---------- new table with      924 representatives

    14329  finished      13778  clusters

Apprixmated maximum memory consumption: 205M
writing new database
writing clustering information
program completed !

Total CPU time 3.83
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. cd-hit.sh). For example:



```

#!/bin/bash
set -e
module load cd-hit
cd-hit -i /fdb/fastadb/drosoph.aa.fas  -o drosoph100 -c 1.00 -n 5 -M 16000 -d 0  -T $SLURM_CPUS_PER_TASK

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] cd-hit.sh
```







