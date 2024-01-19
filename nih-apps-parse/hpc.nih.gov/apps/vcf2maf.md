

document.querySelector('title').textContent = 'vcf2maf on Biowulf';
vcf2maf on Biowulf


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



To convert a VCF into a MAF, each variant must be mapped to only one of all possible gene transcripts/isoforms that it might affect. This selection of a single effect per variant, is often subjective. So this project is an attempt to make the selection criteria smarter, reproducible, and more configurable. And the default criteria must lean towards best practices.



Documentation
* vcf2maf Main Site: [vcf2maf on GitHub](https://github.com/mskcc/vcf2maf)
* Type **vcf2maf.pl --man** to see a brief summary of the options


Important Notes
* Module Name: vcf2maf (see [the modules page](/apps/modules.html) for more information)
 * Multithreaded
 * Environment variables set 
	+ VCF2MAF\_HOME* Example files in $VCF2MAF\_HOME/tests* Reference data in /fdb/VEP/\*/cache -- depends on version of VEP loaded
 * vcf2maf.pl depends heavily on [VEP](VEP.html). At present, only versions of VEP < 107 are compatible.
 * Typically vcf2maf.pl requires at least 20g of memory and minimally 4 cpus.
 * While VEP is tolerant of chromosome format mismatches (when the input .vcf file uses the [UCSC](https://hpc.nih.gov/refdb/dbview.php?id=920) format chr*N* and the reference fasta uses [Ensembl/NCBI](https://hpc.nih.gov/refdb/dbview.php?id=408) format *N*), vcf2maf is not. Make sure the reference fasta chromosome format matches that of your input.



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --cpus-per-task=8 --mem=20g**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load VEP/106 vcf2maf**

[user@cn3144 ~]$ **rm -rf test\_{vcf2maf,maf2maf,maf2vcf}**
[user@cn3144 ~]$ **mkdir test\_vcf2maf && cd $\_**
[user@cn3144 ~]$ **cp $VCF2MAF\_EXAMPLES/test.vcf .**
[user@cn3144 ~]$ **vcf2maf.pl \
 --input-vcf test.vcf \
 --output-maf test.maf \
 --vep-path $VEP\_HOME/bin \
 --vep-data $VEP\_CACHEDIR \
 --ref-fasta $VEP\_CACHEDIR/GRCh38.fa \
 --ncbi-build GRCh38 \
 --cache-version $VEP\_VERSION \
 --vep-forks ${SLURM\_CPUS\_ON\_NODE:-2} \
 --verbose**
[user@cn3144 ~]$ **cd ..**

[user@cn3144 ~]$ **mkdir test\_maf2maf && cd $\_**
[user@cn3144 ~]$ **cp $VCF2MAF\_EXAMPLES/test.maf .**
[user@cn3144 ~]$ **maf2maf.pl \
 --input-maf test.maf \
 --output-maf test.vep.maf \
 --ref-fasta $VEP\_CACHEDIR/GRCh37.fa \
 --vep-path $VEP\_HOME/bin \
 --vep-data $VEP\_CACHEDIR**
[user@cn3144 ~]$ **cd ..**

[user@cn3144 ~]$ **mkdir test\_maf2vcf && cd $\_**
[user@cn3144 ~]$ **cp $VCF2MAF\_EXAMPLES/test.maf .**
[user@cn3144 ~]$ **maf2vcf.pl \
 --input-maf test.maf \
 --output-dir vcfs \
 --ref-fasta $VEP\_CACHEDIR/GRCh37.fa**


```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. vcf2maf.sh). For example:



```
#!/bin/bash
module load VEP/106 vcf2maf
vcf2maf.pl --input-vcf /path/to/data/test.vcf --output-maf test.maf --vep-path $VEP_HOME/bin --vep-data $VEP_CACHEDIR --cache-version $VEP_VERSION \
  --ref-fasta $VEP_CACHEDIR/GRCh38.fa --ncbi-build GRCh38 --vep-forks $SLURM_CPUS_ON_NODE

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] vcf2maf.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. vcf2maf.swarm). For example:



```
vcf2maf.pl --input-vcf data/test1.vcf --output-maf data/test1.maf --vep-path $VEP_HOME/bin --vep-data $VEP_CACHEDIR --cache-version $VEP_VERSION --ref-fasta $VEP_CACHEDIR/GRCh38.fa --ncbi-build GRCh38 --vep-forks $SLURM_CPUS_ON_NODE
vcf2maf.pl --input-vcf data/test2.vcf --output-maf data/test2.maf --vep-path $VEP_HOME/bin --vep-data $VEP_CACHEDIR --cache-version $VEP_VERSION --ref-fasta $VEP_CACHEDIR/GRCh38.fa --ncbi-build GRCh38 --vep-forks $SLURM_CPUS_ON_NODE
vcf2maf.pl --input-vcf data/test3.vcf --output-maf data/test3.maf --vep-path $VEP_HOME/bin --vep-data $VEP_CACHEDIR --cache-version $VEP_VERSION --ref-fasta $VEP_CACHEDIR/GRCh38.fa --ncbi-build GRCh38 --vep-forks $SLURM_CPUS_ON_NODE
vcf2maf.pl --input-vcf data/test4.vcf --output-maf data/test4.maf --vep-path $VEP_HOME/bin --vep-data $VEP_CACHEDIR --cache-version $VEP_VERSION --ref-fasta $VEP_CACHEDIR/GRCh38.fa --ncbi-build GRCh38 --vep-forks $SLURM_CPUS_ON_NODE
```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f vcf2maf.swarm [-g #] [-t #] --module vcf2maf VEP/106
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module vcf2maf Loads the vcf2maf module for each subjob in the swarm 
 | |
 | |
 | |








