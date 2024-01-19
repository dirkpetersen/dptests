

document.querySelector('title').textContent = 'FAN-C on Biowulf';
FAN-C on Biowulf


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



FAN-C provides a pipeline for analysing Hi-C data starting at mapped paired-end sequencing reads.



### Reference:


* [Kruse, Kai, Clemens B. Hug, and Juan M. Vaquerizas. "FAN-C: a feature-rich framework for the analysis and visualisation of chromosome conformation capture data." *Genome biology* 21.1 (2020): 1-19.](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02215-9)


Documentation
* [FAN-C Main Site](https://fan-c.readthedocs.io/en/latest/)
* [FAN-C on GitHub](https://github.com/vaquerizaslab/fanc)


Important Notes
* Module Name: fanc (see [the modules page](/apps/modules.html) for more information)
 * FAN-C is multithreaded. Please use the option argument pair -t $SLURM\_CPUS\_PER\_TASK as shown below. 
 * This application produces PDF reports. You can use [hpcdrive](/docs/hpcdrive.html) to view these reports on your local workstation.
* Path to example files is set in the FANC\_EXAMPLES environment variable after loading the module.



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf ~]$ **sinteractive -c8 --mem=16g --gres=lscratch:30**
salloc.exe: Pending job allocation 11754716
salloc.exe: job 11754716 queued and waiting for resources
salloc.exe: job 11754716 has been allocated resources
salloc.exe: Granted job allocation 11754716
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn0996 are ready for job
srun: error: x11: no local DISPLAY defined, skipping
error: unable to open file /tmp/slurm-spank-x11.11754716.0
slurmstepd: error: x11: unable to read DISPLAY value

[user@cn0996 ~]$ **module load fanc**
[+] Loading fanc  0.9.17  on cn0996
[+] Loading singularity  3.7.2  on cn0996

[user@cn0996 ~]$ **cd /lscratch/$SLURM\_JOB\_ID**

[user@cn0996 11754716]$ **cp $FANC\_EXAMPLES .**

[user@cn0996 11754716]$ **unzip examples.zip**
Archive:  examples.zip
   creating: examples/
  inflating: examples/SRR4271982_chr18_19_2.fastq.gzip
[snip...]
  inflating: examples/SRR4271982_chr18_19_1.fastq.gzip

[user@cn0996 11754716]$ **cd examples/**

[user@cn0996 examples]$ **fanc auto SRR4271982\_chr18\_19\_1.fastq.gzip SRR4271982\_chr18\_19\_2.fastq.gzip \
 output/ -g hg19\_chr18\_19.fa -i hg19\_chr18\_19/hg19\_chr18\_19 -n fanc\_example \
 -t ${SLURM\_CPUS\_PER\_TASK} -r HindIII --split-ligation-junction -q 30**
2021-03-30 13:07:51,440 INFO FAN-C version: 0.9.17
2021-03-30 13:07:51,443 INFO Output folder: output/
2021-03-30 13:07:51,443 INFO Input files: SRR4271982_chr18_19_1.fastq.gzip, SRR4271982_chr18_19_2.fastq.gzip
2021-03-30 13:07:51,443 INFO Input file types: fastq, fastq
2021-03-30 13:07:51,443 INFO Final basename: fanc_example (you can change this with the -n option!)
[...snip]
2021-03-30 13:22:16,001 INFO Total: 556663. Filtered: 294
Expected 100% (556369 of 556369) |###########################################| Elapsed Time: 0:00:08 Time:  0:00:08
2021-03-30 13:22:29,543 INFO Done.
2021-03-30 13:22:33,641 INFO Saving statistics...
2021-03-30 13:22:33,832 INFO Normalising binned Hic file
Expected 100% (556369 of 556369) |###########################################| Elapsed Time: 0:00:08 Time:  0:00:08
Closing remaining open files:output/hic/fanc_example.hic...done

[user@cn0996 examples]$ **fancplot chr18:63mb-70mb -p triangular -vmax 0.05 output/hic/binned/fanc\_example\_100kb.hic**
2021-03-30 13:29:15,888 INFO Found 1 regions
Closing remaining open files:output/hic/binned/fanc_example_100kb.hic...done

[user@cn0996 examples]$ **ls -lh**
total 558M
drwxr-xr-x 9 user user 4.0K Aug  6  2020 architecture
drwxr-xr-x 2 user user 4.0K Feb  4  2020 bwa-index
drwxr-xr-x 2 user user 4.0K Feb  4  2020 hg19_chr18_19
-rw-r--r-- 1 user user 134M Feb  4  2020 hg19_chr18_19.fa
-rw-r--r-- 1 user user 996K Feb  4  2020 hg19_chr18_19_re_fragments.bed
drwxr-xr-x 2 user user 4.0K Aug  6  2020 hicpro
drwxr-x--- 7 user user 4.0K Mar 30 13:07 output
-rw-r--r-- 1 user user 200M Feb  4  2020 SRR4271982_chr18_19_1.fastq.gzip
-rw-r--r-- 1 user user 200M Feb  4  2020 SRR4271982_chr18_19_2.fastq.gzip
-rw-r--r-- 1 user user  52K Feb  4  2020 test.cool
-rw-r--r-- 1 user user  24M Feb  4  2020 test.hic

[user@cn0996 examples]$ **exit**
exit
salloc.exe: Relinquishing job allocation 11754716
salloc.exe: Job allocation 11754716 has been revoked.

[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. fanc.sh). For example:



```

#!/bin/bash
set -e
module load fanc

fanc auto SRR4271982_chr18_19_1.fastq.gzip SRR4271982_chr18_19_2.fastq.gzip \
    output/ -g hg19_chr18_19.fa -i hg19_chr18_19/hg19_chr18_19 -n fanc_example \
    -t ${SLURM_CPUS_PER_TASK} -r HindIII --split-ligation-junction -q 30

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] fanc.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. fanc.swarm). For example:



```

fanc auto a.fastq.gzip chrN.fastq.gzip outa -g N.fa -t $SLURM_CPUS_PER_TASK <opts> 
fanc auto b.fastq.gzip chrN.fastq.gzip outb -g N.fa -t $SLURM_CPUS_PER_TASK <opts> 
fanc auto c.fastq.gzip chrN.fastq.gzip outc -g N.fa -t $SLURM_CPUS_PER_TASK <opts> 
fanc auto c.fastq.gzip chrN.fastq.gzip outd -g N.fa -t $SLURM_CPUS_PER_TASK <opts> 

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f fanc.swarm [-g #] [-t #] --module fanc
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module fanc Loads the fanc module for each subjob in the swarm 
 | |
 | |
 | |








