

document.querySelector('title').textContent = 'Vcftools on Biowulf';
Vcftools on Biowulf


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


VCFtools contains a Perl API (Vcf.pm) and a number of Perl scripts that 
 can be used to perform common tasks with VCF files such as file validation, 
 file merging, intersecting, complements, etc. The Perl tools support all 
 versions of the VCF specification (3.2, 3.3, and 4.0), nevertheless, the 
 users are encouraged to use the latest version VCFv4.0. The VCFtools in 
 general have been used mainly with diploid data, but the Perl tools aim 
 to support polyploid data as well.


VCFTools is maintained and developed by Adam Auton, Peter Danecek and collaborators.


  




### References:


* vcftools paper


Documentation
* <https://vcftools.github.io/examples.html>


Important Notes * Module Name: vcftools (see [the 
 modules page](/apps/modules.html) for more information)



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load vcftools**
[user@cn3144 ~]$ **compare-vcf inputFile1 inputFile2**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. vcftools.sh). For example:



```

#!/bin/bash
set -e
module load vcftools
compare-vcf inputFile1 inputFile2
```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --mem=10g vcftools.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. vcftools.swarm). For example:



```
cd dir1;compare-vcf inputFile1 inputFile2
cd dir2;compare-vcf inputFile1 inputFile2
cd dir3;compare-vcf inputFile1 inputFile2
cd dir4;compare-vcf inputFile1 inputFile2
cd dir5;compare-vcf inputFile1 inputFile2
```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f vcftools.swarm -g 10 --module vcftools
```

where
 

|  |  |  |  |
| --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in 
 the swarm command file) 
 | --module vcftools Loads the fastqc module for each subjob in the swarm 
  | |
 | |








