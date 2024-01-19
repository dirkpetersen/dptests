

document.querySelector('title').textContent = 'muTect on Biowulf';
muTect on Biowulf


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



MuTect is a method developed at the Broad Institute for the reliable and accurate identification of somatic point mutations in next generation sequencing data of cancer genomes.



### References:




Cibulskis, K. et al.
[**Sensitive detection of somatic point mutations in impure and heterogeneous cancer samples.**](https://www.ncbi.nlm.nih.gov/pubmed/23396013)
*Nat Biotechnology (2013) Mar; 31(3)213-9* * Paper


**NOTE:** MuTect has been merged into the most recent versions of GATK, and is no longer 
supported as a separate application. Please see [https://hpc.nih.gov/apps/GATK.html](GATK.html) for
information on using GATK.


Documentation
* [MuTect Homepage](http://www.broadinstitute.org/cancer/cga/mutect/)
* [MuTect Discussion Forum](http://gatkforums.broadinstitute.org/categories/mutect)


Important Notes
* Module Name: muTect (see [the modules page](/apps/modules.html) for more information)
* Singlethreaded
* Environment variables set 
	+ MUTECT\_JAR - full path to the mutect jar
	+ MUTECT\_JARPATH - path to the directory containing the mutect jar* Reference data in /fdb/muTect/



**NOTE:** muTect uses code base from GATK, and therefore has many of the same options. One option,
 **-nt** or **--num\_threads** **DOES NOT** work properly. **DO NOT** use this
 option.


MuTect requires two BAM input files, one for normal tissues, the other for the tumor tissue. MuTect outputs
 a wiggle format coverage file. An additional wiggle file can be generated to display observed depth.


Two extra options have been added to allow for memory allocation and temporary file directory.


* **--memory** memory allocated (default = 2gb)
* **--tmpdir** tmpdir location (default = /lscratch)


By default, muTect uses 2gb of memory. To allocate
 5gb of memory, include **--memory 5g** on the commandline.


MuTect takes as parameters database files, depending on the build of your alignments and which dbSNP version
 you are using. These files are located in **/fdb/muTect**.


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

[user@cn3144 ~]$

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. muTect.sh). For example:



```

#!/bin/bash
module load muTect
muTect --memory 8g --analysis_type MuTect \
--reference_sequence /fdb/muTect/ucsc.hg19.fasta \
--dbsnp              /fdb/muTect/dbsnp_137.hg19.vcf \
--cosmic             /fdb/muTect/cosmic_v67.hg19.vcf \
--input_file:normal /full/path/to/Normal.cleaned.bam \
--input_file:tumor /full/path/to/Tumor.cleaned.bam \
--out /full/path/to/example.call_stats.txt \ 
--coverage_file /full/path/to/example.coverage.wig.txt

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] muTect.sh
```







