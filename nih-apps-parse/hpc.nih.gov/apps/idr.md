

document.querySelector('title').textContent = 'IDR on Biowulf';
IDR on Biowulf


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



The IDR (Irreproducible Discovery Rate) framework is a uniﬁed approach to measure the reproducibility of ﬁndings identiﬁed from replicate experiments and provide highly stable thresholds based on reproducibility. Unlike the usual scalar measures of reproducibility, the IDR approach creates a curve, which quantitatively assesses when the ﬁndings are no longer consistent across replicates. In layman's terms, the IDR method compares a pair of ranked lists of identifications (such as ChIP-seq peaks). These ranked lists should not be pre-thresholded i.e. they should provide identifications across the entire spectrum of high confidence/enrichment (signal) and low confidence/enrichment (noise). The IDR method then fits the bivariate rank distributions over the replicates in order to separate signal from noise based on a defined confidence of rank consistency and reproducibility of identifications i.e the IDR threshold.



### References:


* "Measuring reproducibility of high-throughput experiments" (2011), Annals of Applied Statistics, Vol. 5, No. 3, 1752-1779, by Li, Brown, Huang, and Bickel


Documentation
* [IDR Main Site](https://github.com/nboley/idr)


Important Notes
* Module Name: idr (see [the modules page](/apps/modules.html) for more information)
* Example files in /data/idr



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load idr**
[+] Loading idr, version 2.0.3... 
[user@cn3144 ~]$ **idr --samples /fdb/idr/peak{1,2}**
Initial parameter values: [0.10 1.00 0.20 0.50]
Final parameter values: [1.57 1.26 0.89 0.41]
Number of reported peaks - 50537/50537 (100.0%)

Number of peaks passing IDR cutoff of 0.05 - 12748/50537 (25.2%)

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. idr.bash). For example:



```

#!/bin/bash
set -e
module load idr
idr --samples /fdb/idr/peak{1,2} -o peak1-2.txt

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] idr.bash
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. idr.swarm). For example:



```

idr --samples set1/peak{1,2} -o idr/set1.txt
idr --samples set2/peak{1,2} -o idr/set2.txt
idr --samples set3/peak{1,2} -o idr/set3.txt
idr --samples set4/peak{1,2} -o idr/set4.txt

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f idr.swarm [-g #] [-t #] --module idr
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module idr Loads the idr module for each subjob in the swarm 
 | |
 | |
 | |








