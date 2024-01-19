

document.querySelector('title').textContent = "slivar";
slivar on Biowulf


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



 slivar is a set of command-line tools that enables rapid querying and filtering of VCF files. It facilitates operations on trios and groups and allows arbitrary expressions using simple javascript.




 use cases for slivar:
 * annotate variants with gnomad allele frequencies from combined exomes + whole genomes at > 30K variants/second using only a 1.5GB compressed annotation file.
* call denovo variants with a simple expression that uses mom, dad, kid labels that is applied to each trio in a cohort (as inferred from a pedigree file). kid.het && mom.hom\_ref && dad.hom\_ref && kid.DP > 10 && mom.DP > 10 && dad.DP > 10
* define and filter on arbitrary groups with labels. For example, 7 sets of samples each with 1 normal and 3 tumor time-points: normal.AD[0] = 0 && tumor1.AB < tumor2.AB && tumor2.AB < tumor3.AB
* filter variants with simple expressions: variant.call\_rate > 0.9 && variant.FILTER == "PASS" && INFO.AC < 22 && variant.num\_hom\_alt == 0
* see [using slivar for rare disease research](https://github.com/brentp/slivar/wiki/rare-disease)





### References:


* Pedersen, B.S., Brown, J.M., Dashnow, H. et al.
 [**Effective variant filtering and expected candidate variant yield in studies of rare human disease.**](https://www.ncbi.nlm.nih.gov/pubmed/00000000)
 *npj Genom. Med. 6, 60 (2021).*
[doi:10.1038/s41525-021-00227-3](https://doi.org/10.1038/s41525-021-00227-3)


Documentation
* [slivar wiki](https://github.com/brentp/slivar/wiki)
* [slivar on GitHub](https://github.com/brentp/slivar)


Important Notes
* Module Name: slivar (see [the modules page](/apps/modules.html) for more information)
 * For multithreaded application of the expr subcommand only, use [pslivar](https://github.com/brentp/slivar/wiki/parallel-slivar).
* Environment variables set 
	+ SLIVAR\_HOME* Developer's provided javascript expressions in $SLIVAR\_HOME/js
* Example files in $SLIVAR\_HOME/tests
* Reference data (pregenerated gnotation files) in /fdb/slivar/



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

[user@cn3144 ~]$ **module load slivar**

[user@cn3144 ~]$ **slivar expr \
 --js $SLIVAR\_HOME/js/slivar-functions.js \
 -g /fdb/slivar/gnomad.hg38.genomes.v3.fix.zip \
 --vcf $SLIVAR\_HOME/tests/ashk-trio.vcf.gz \
 --ped $SLIVAR\_HOME/tests/ashk-trio.ped \
 --info "INFO.gnomad\_popmax\_af < 0.01 && variant.FILTER == 'PASS'" \
 --trio "example\_denovo:denovo(kid, dad, mom)" \
 --family-expr "denovo:fam.every(segregating\_denovo)" \
 --trio "custom:kid.het && mom.het && dad.het && kid.GQ > 20 && mom.GQ > 20 && dad.GQ > 20" \
 --pass-only \
 -o ashk-trio.slivar.vcf**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. slivar.sh). For example:



```

#!/bin/bash
set -e
module load slivar
slivar expr \
  --js $SLIVAR_HOME/js/slivar-functions.js \
  -g /fdb/slivar/gnomad.hg38.genomes.v3.fix.zip \
  --vcf $SLIVAR_HOME/tests/ashk-trio.vcf.gz \
  --ped $SLIVAR_HOME/tests/ashk-trio.ped \
  --info "INFO.gnomad_popmax_af < 0.01 && variant.FILTER == 'PASS'" \
  --trio "example_denovo:denovo(kid, dad, mom)" \
  --family-expr "denovo:fam.every(segregating_denovo)" \
  --trio "custom:kid.het && mom.het && dad.het && kid.GQ > 20 && mom.GQ > 20 && dad.GQ > 20"  \
  --pass-only \
  -o ashk-trio.slivar.vcf

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] slivar.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. slivar.swarm). For example:



```

slivar expr --js $SLIVAR_HOME/js/slivar-functions.js -g /fdb/slivar/gnomad.hg38.genomes.v3.fix.zip --vcf sample1.vcf --ped sample1.ped -o sample1.slivar.vcf ...
slivar expr --js $SLIVAR_HOME/js/slivar-functions.js -g /fdb/slivar/gnomad.hg38.genomes.v3.fix.zip --vcf sample2.vcf --ped sample2.ped -o sample2.slivar.vcf ...
slivar expr --js $SLIVAR_HOME/js/slivar-functions.js -g /fdb/slivar/gnomad.hg38.genomes.v3.fix.zip --vcf sample3.vcf --ped sample3.ped -o sample3.slivar.vcf ...
slivar expr --js $SLIVAR_HOME/js/slivar-functions.js -g /fdb/slivar/gnomad.hg38.genomes.v3.fix.zip --vcf sample4.vcf --ped sample4.ped -o sample4.slivar.vcf ...

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f slivar.swarm [-g #] [-t #] --module slivar
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module slivar Loads the slivar module for each subjob in the swarm 
 | |
 | |
 | |








