

document.querySelector('title').textContent = 'Ultraplex on Biowulf';
Ultraplex on Biowulf


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



Ultraplex is an application for demultiplexing and processing FASTQ files. The processing steps include removing low quality bases, removing sequencing adaptors, and separating FASTQ using demultiplexing barcodes. Ultraplex is intended for use with custom library preparation protocols instead of commercial prep kits.



Documentation
* [Ultraplex Github](https://github.com/ulelab/ultraplex)


Important Notes
This application requires a [lscratch allocation](/docs/userguide.html#local)


* Module Name: ultraplex (see [the modules page](/apps/modules.html) for more information)
* Multithreaded
* Environment variables set 
	+ ULTRAPLEX\_TEST\_DATA* Test data files in /usr/local/apps/ultraplex/TEST\_DATA



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive --gres=lscratch:10**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load ultraplex**
[+] Loading singularity  3.8.5-1  on cn0847
[+] Loading ultraplex 1.2.5  ...

[user@cn3144 ~]$ **cd /lscratch/$SLURM\_JOB\_ID**

[user@cn3144 ~]$ **cp -r $ULTRAPLEX\_TEST\_DATA tests**

[user@cn3144 ~]$ **ultraplex -i tests/reads1.fastq.gz \
 -i2 tests/reads2.fastq.gz \
 -b tests/barcodes\_5\_and\_3.csv \
 -d PE -o paired\_end**
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
@@@@@   @@@@   .@@   @@@@@@@          @@         @@@@@@    (@@@@@        @@@    @@@@@@@         @@    @@@    @
@@@@   @@@@    @@   ,@@@@@@@@@    @@@@@    @@@   @@@@   (   @@@@   @@@   @@    @@@@@@@   @@@@@@@@@@   #   @@@@
@@@   &@@@    @@    @@@@@@@@@%   @@@@@         @@@@(   @    @@@         @@    @@@@@@@         @@@@@     @@@@@@
@@    @@@    @@    @@@@@@@@@@   @@@@@    @@    @@@          @@   .@@@@@@@*   @@@@@@@    @@@@@@@@.   @    @@@@@
@@@       @@@@.        @@@@@   @@@@@&   @@@.   &   &@@@@    @    @@@@@@@@         @          @    @@@@    @@@@
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

Namespace(adapter='AGATCGGAAGAGCGGTTCAG', adapter2='AGATCGGAAGAGCGTCGTG', barcodes='tests/barcodes_5_and_3.csv', directory='PE', dont_build_reference=False, final_min_length=20, fiveprimemismatches=1, ignore_no_match=False, ignore_space_warning=False, input_2='tests/reads2.fastq.gz', inputfastq='tests/reads1.fastq.gz', keep_barcode=False, min_trim=3, outputprefix='paired_end', phredquality=30, phredquality_5_prime=0, sbatchcompression=False, threads=4, threeprimemismatches=0, ultra=False)
Demultiplexing...
Demultiplexing complete! 2700 reads processed in 0.0 seconds

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. ultraplex.sh). For example:



```

#!/bin/bash
set -e
module load ultraplex
ultraplex -i R1.fastq.gz \
          -i2 R2.fastq.gz \
          -b barcods.csv \
          -d PE \
          -t $SLURM_CPUS_PER_TASK
          -o output_dir

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] --gres=lscratch:# ultraplex.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. ultraplex.swarm). For example:



```

ultraplex -t $SLURM_CPUS_PER_TASK -i L1.R1.fastq.gz -i2 L1.R2.fastq.gz -b barcodes.csv -d PE -o L1
ultraplex -t $SLURM_CPUS_PER_TASK -i L2.R1.fastq.gz -i2 L2.R2.fastq.gz -b barcodes.csv -d PE -o L2
ultraplex -t $SLURM_CPUS_PER_TASK -i L3.R1.fastq.gz -i2 L3.R2.fastq.gz -b barcodes.csv -d PE -o L3
ultraplex -t $SLURM_CPUS_PER_TASK -i L4.R1.fastq.gz -i2 L4.R2.fastq.gz -b barcodes.csv -d PE -o L4

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f ultraplex.swarm [-g #] [-t #] --gres=lscratch:# --module ultraplex
```

where


|  |  |
| --- | --- |
| -g *#* | Number of Gigabytes of memory required for each process (1 line in the swarm command file) |
| -t *#* | Number of threads/CPUs required for each process (1 line in the swarm command file). |
| --gres=lscratch:# | lscratch amount in GB allocated for each process (1 line in the swarm command file). |
| --module ultraplex | Loads the Ultraplex module for each subjob in the swarm |








