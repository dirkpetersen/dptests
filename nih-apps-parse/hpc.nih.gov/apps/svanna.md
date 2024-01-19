
















svanna on Biowulf


  









|  |  |  |
| --- | --- | --- |
|  | 
[Biowulf High Performance Computing at the NIH](https://hpc.nih.gov)
 | 








[GitHub](https://github.com/NIH-HPC)
[YouTube](https://www.youtube.com/channel/UCx-kNd1kBskYr5KLT9-Erew)
[@nih_hpc](https://twitter.com/nih_hpc)
[RSS Feed](/hpc_RSS.xml)
 |




* [Systems](https://hpc.nih.gov/systems/)
* [Applications](https://hpc.nih.gov/apps/)
* [Reference Data](https://hpc.nih.gov/refdb/)
* [Storage](https://hpc.nih.gov/storage/)
* [User Guides](https://hpc.nih.gov/docs/user_guides.html)
* [Training](https://hpc.nih.gov/training/)
* [User Dashboard](https://hpcnihapps.cit.nih.gov/auth/dashboard/)
* [How To](https://hpc.nih.gov/docs/how_to.html)
* [About](https://hpc.nih.gov/about/)







svanna on Biowulf


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



The svanna is an efficient and accurate pathogenicity prediction for coding and regulatory structural variants in long-read genome sequencing.






### References:


* Danis, D., Jacobsen, J.O.B., Balachandran, P. et al. *SvAnna: efficient and accurate pathogenicity prediction of coding and regulatory structural variants in long-read genome sequencing.*Genome Med 14, 44 (2022)

 [Journal](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-022-01046-6)


Documentation
* svanna Github:[Github](https://github.com/TheJacksonLaboratory/SvAnna)
* svanna doc:[Github](https://svanna.readthedocs.io/en/master/)


Important Notes
* Module Name: svanna (see [the modules page](/apps/modules.html) for more information)
 * The svanna ((available on biowulf) command lines could be run as:
 
```

	java -jar $SVANNA_PATH/svanna-cli-1.0.3.jar prioritize -d $SVANNA_DB/2304 --help
	
```
* svanna database in $SVANNA\_DB 
 * svanna testing vcf in $SVANNA\_TEST\_DATA


Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]8$ **sinteractive --mem=8g**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load svanna**
[user@cn3144 ~]$ **mkdir /data/$USER/svanna\_test/**
[user@cn3144 ~]$ **cd /data/$USER/svanna\_test/**
[user@cn3144 ~]$ **java -jar $SVANNA\_PATH/svanna-cli-1.0.3.jar prioritize -d $SVANNA\_DB/2304 --help**
Prioritize the variants.
Usage: svanna-cli.jar prioritize [-hVv] [--no-breakends] [--uncompressed-output] -d=path/to/datadir
                                 [--frequency-threshold=] [--ic-mica-mode={DATABASE,IN\_MEMORY}]
 [--min-read-support=] [--n-threads=2] [--out-dir=]
 [--output-format=html] [--overlap-threshold=] [-p=]
 [--prefix=] [--promoter-fitness-gain=]
 [--promoter-length=] [--report-top-variants=100]
 [--term-similarity-measure={RESNIK\_SYMMETRIC, RESNIK\_ASYMETRIC}] [--vcf=]
 [-t=]...
 -v Specify multiple -v options to increase verbosity.
 For example, `-v -v -v` or `-vvv`
 -d, --data-directory=path/to/datadir
 Path to SvAnna data directory.
 -h, --help Show this help message and exit.
 -V, --version Print version information and exit.
SvAnna configuration:
 --term-similarity-measure={RESNIK\_SYMMETRIC, RESNIK\_ASYMETRIC}
 Phenotype term similarity measure (default: RESNIK\_SYMMETRIC).
 --ic-mica-mode={DATABASE,IN\_MEMORY}
 The mode for getting information content of the most informative common ancestors for
 terms t1, and t2 (default: DATABASE).
 --promoter-length=
 Number of bases prepended to a transcript and evaluated as a promoter region (default:
 2000).
 --promoter-fitness-gain=
 Set to 0. to score the promoter variants as strictly as coding variants, or to 1. to skip
 them altogether (default: 0.6).
Analysis input:
 -p, --phenopacket=
 Path to v1 or v2 phenopacket in JSON, YAML or Protobuf format.
 -t, --phenotype-term=
 HPO term ID(s). Can be provided multiple times.
 --vcf= Path to the input VCF file.
Run options:
 --frequency-threshold=
 Frequency threshold as a percentage [0-100] (default: 1.0).
 --overlap-threshold=
 Percentage threshold for determining variant's region is similar enough to database entry
 (default: 80.0).
 --min-read-support=
 Minimum number of ALT reads to prioritize (default: 3).
 --n-threads=2 Process variants using n threads (default: 2).
Output options:
 --no-breakends Do not include breakend variants into HTML report (default: false).
 --output-format=html Comma separated list of output formats to use for writing the results (default: html).
 --out-dir= Path to folder where to write the output files (default: current working directory).
 --prefix= Prefix for output files (default: based on the input VCF name).
 --report-top-variants=100
 Report top n variants (default: 100).
 --uncompressed-output Write tabular and VCF output formats with no compression (default: false).
See the full documentation at `https://svanna.readthedocs.io/en/master`

[user@cn3144 ]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. svanna.sh). For example:




hljs.highlightAll();

```


#!/bin/bash
set -e
module load svanna
cp $SVANNA_TEST_DATA/*.vcf .
java -jar $SVANNA_PATH/svanna-cli-1.0.3.jar prioritize -d $SVANNA_DB/2304 --vcf example.vcf  --phenotype-term HP:0011890 --phenotype-term HP:0000978 --phenotype-term HP:0012147


```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --mem=8g svanna.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. svanna.swarm). For example:



```

java -jar $SVANNA_PATH/svanna-cli-1.0.3.jar prioritize -d $SVANNA_DB/2304 --vcf example1.vcf --out-dir results --prefix example1  --phenotype-term HP:0011890 --phenotype-term HP:0000978 --phenotype-term HP:0012147
java -jar $SVANNA_PATH/svanna-cli-1.0.3.jar prioritize -d $SVANNA_DB/2304 --vcf example2.vcf --out-dir results --prefix example2  --phenotype-term HP:0011890 --phenotype-term HP:0000978 --phenotype-term HP:0012147

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f svanna.swarm [-t #] [-g #] --module svanna
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module svanna Loads the svanna module for each subjob in the swarm
 | |
 | |
 | |











[HPC @ NIH](https://hpc.nih.gov)  ~
[Contact](https://hpc.nih.gov/about/contact.html)


[Disclaimer](https://hpc.nih.gov/docs/disclaimer.html) ~ 
[Privacy](https://hpc.nih.gov/docs/privacy.html) ~ 
[Accessibility](https://hpc.nih.gov/docs/accessibility.html) ~ 
[CIT](https://cit.nih.gov/) ~ 
[NIH](https://www.nih.gov/) ~ 
[DHHS](https://www.dhhs.gov/) ~ 
[USA.gov](https://www.firstgov.gov/) ~
[HHS Vulnerability Disclosure](https://www.hhs.gov/vulnerability-disclosure-policy/index.html)



  


