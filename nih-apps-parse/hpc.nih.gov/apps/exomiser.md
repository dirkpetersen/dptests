

document.querySelector('title').textContent = 'Exomiser on Biowulf';
Exomiser on Biowulf


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


The Exomiser is a Java program that finds potential disease-causing variants 
 from whole-exome or whole-genome sequencing data.


Starting from a VCF file and a set of phenotypes encoded using the [Human 
 Phenotype Ontology (HPO)](http://www.human-phenotype-ontology.org/) it will annotate, filter and prioritise likely 
 causative variants. The program does this based on user-defined criteria 
 such as a variant's predicted pathogenicity, frequency of occurrence in 
 a population and also how closely the given phenotype matches the known 
 phenotype of diseased genes from human and model organism data.


The functional annotation of variants is handled by [Jannovar](https://github.com/charite/jannovar) 
 and uses [UCSC](http://genome.ucsc.edu/) KnownGene transcript 
 definitions and hg19 genomic coordinates.


Variants are prioritised according to user-defined criteria on variant 
 frequency, pathogenicity, quality, inheritance pattern, and model organism 
 phenotype data. Predicted pathogenicity data is extracted from the [dbNSFP](http://www.ncbi.nlm.nih.gov/pubmed/21520341) 
 resource. Variant frequency data is taken from the [1000 
 Genomes](http://www.1000genomes.org/), [ESP](http://evs.gs.washington.edu/EVS) and [ExAC](http://exac.broadinstitute.org/) 
 datasets. Subsets of these frequency and pathogenicity data can be defined 
 to further tune the analysis. Cross-species phenotype comparisons come from 
 our PhenoDigm tool powered by the OWLTools [OWLSim](https://github.com/owlcollab/owltools) 
 algorithm.


### References:


* The Exomiser was developed by the Computational Biology and Bioinformatics 
 group at the Institute for Medical Genetics and Human Genetics of the 
 Charité - Universitätsmedizin Berlin, the Mouse Informatics Group at the 
 Sanger Institute and other members of the [Monarch 
 initiative](https://monarchinitiative.org/).


Documentation
* <https://github.com/exomiser/Exomiser>


Important Notes
* Module Name: exomiser (see [the 
 modules page](/apps/modules.html) for more information) 
 * Example files can be copied from /usr/local/apps/exomiser/examples/
 * applications.properties must exist in current working directory 
 * exomiser tends to auto thread to ~33 threads. So --cpus-per-task=34 is recommended
 * output directory has to exist for the job to proceed or job may terminated with commandlinerunner error
 * For a first run, you will likely have no idea how much mmemory is required. Submit a test job requesting, say, 20 GB of memory, and set the java process to use, say, 15 or 18 GB of memory (-Xms2g -Xmx18g). Xms indicates the initial memory allocation, and Xmx indicates the maximum memory allocation. If your job completes successfully, then look at the dashboard or use jobhist to see how much memory was actually used, and use that value (with a 10% buffer) to submit future jobs.



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. 
In the example below, previous testing has shown that the exomiser process on this test data requires about 3.3 GB to run. Therefore the job is submitted requesting 5 GB of memory, and the java process requests 4 GB (-Xmx4g in the command below). For other test data, more memory may be required. 



Sample session: (user input in bold)



```

[user@biowulf]$ **sinteractive --mem=5g**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load exomiser**
[user@cn3144 ~]$ **mkdir /data/$USER/exomiser**[user@cn3144 ~]$ **cd /data/$USER/exomiser**[user@cn3144 ~]$ **cp -rp /usr/local/apps/exomiser/9.0.0/examples .**[user@cn3144 ~]$ **cp /usr/local/apps/exomiser/9.0.0/application.properties .**[user@cn3144 ~]$ **mkdir results**[user@cn3144 ~]$ **java -Xms2g -Xmx4g -jar $EXOMISERJAR --analysis examples/test-analysis-exome.yml**

08:57:33.907 [main] INFO org.monarchinitiative.exomiser.cli.Main - Locale set to en_GB

 Welcome to:
  _____ _            _____                     _
 |_   _| |__   ___  | ____|_  _____  _ __ ___ (_)___  ___ _ __
   | | | '_ \ / _ \ |  _| \ \/ / _ \| '_ ` _ \| / __|/ _ \ '__|
   | | | | | |  __/ | |___ >  < (_) | | | | | | \__ \  __/ |
   |_| |_| |_|\___| |_____/_/\_\___/|_| |_| |_|_|___/\___|_|

 A Tool to Annotate and Prioritize Exome Variants     v11.0.0

2019-08-09 08:57:35.487  INFO 14526 --- [           main] org.monarchinitiative.exomiser.cli.Main  : Starting Main on cn3096 with PID 14526 (/usr/local/apps/exomiser/11.0.0/exomiser-cli-11.0.0.jar started by user in /spin1/users/user/exomiser)
2019-08-09 08:57:35.496  INFO 14526 --- [           main] org.monarchinitiative.exomiser.cli.Main  : No active profile set, falling back to default profiles: default
2019-08-09 08:57:35.606  INFO 14526 --- [           main] s.c.a.AnnotationConfigApplicationContext : Refreshing org.springframework.context.annotation.AnnotationConfigApplicationContext@47542153: startup date [Fri Aug 09 08:57:35 EDT 2019]; root of context hierarchy
2019-08-09 08:57:40.127  INFO 14526 --- [           main] o.m.exomiser.cli.config.MainConfig       : Exomiser home: /usr/local/apps/exomiser/11.0.0
2019-08-09 08:57:40.222  INFO 14526 --- [           main] o.m.exomiser.cli.config.MainConfig       : Data source directory defined in properties as: /usr/local/apps/exomiser/9.0.0/data
2019-08-09 08:57:40.250  INFO 14526 --- [           main] o.m.exomiser.cli.config.MainConfig       : Root data source directory set to: /usr/local/apps/exomiser/9.0.0/data
2019-08-09 08:57:40.383  INFO 14526 --- [           main] com.zaxxer.hikari.HikariDataSource       : exomiser-phenotype-1711 - Starting...
2019-08-09 08:57:42.380  INFO 14526 --- [           main] com.zaxxer.hikari.HikariDataSource       : exomiser-phenotype-1711 - Start completed.
2019-08-09 08:57:50.294  INFO 14526 --- [           main] o.m.e.c.prioritisers.util.DataMatrixIO   : reading line 500
[...]

2019-08-09 09:01:07.738  INFO 14526 --- [           main] com.zaxxer.hikari.HikariDataSource       : exomiser-phenotype-1711 - Shutdown initiated...
2019-08-09 09:01:07.775  INFO 14526 --- [           main] com.zaxxer.hikari.HikariDataSource       : exomiser-phenotype-1711 - Shutdown completed.
2019-08-09 09:01:07.776  INFO 14526 --- [           main] org.monarchinitiative.exomiser.cli.Main  : Exomising finished - Bye!

[user@cn3144 ~]$ **ls -rtl results**
total 6476
-rw-r--r-- 1 user user 5764102 Aug  9 09:01 Pfeiffer-hiphive-exome-SPARSE.html
-rw-r--r-- 1 user user  444683 Aug  9 09:01 Pfeiffer-hiphive-exome-SPARSE.vcf
-rw-r--r-- 1 user user  161878 Aug  9 09:01 Pfeiffer-hiphive-exome-SPARSE.genes.tsv
-rw-r--r-- 1 user user  212547 Aug  9 09:01 Pfeiffer-hiphive-exome-SPARSE.variants.tsv

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. exomiser.sh). For example:



```

#!/bin/bash
set -e
module load exomiser
cd /data/$USER/..........

java -Xms2g -Xmx10g -jar $EXOMISERJAR --analysis examples/test-analysis-exome.yml
```

NOTE: -Xmx10g is depend on how much memory is requested when submitting 
 the job below. For example, if 10g of memory is requested, then put 10g 
 here. 


Submit this job using the Slurm [sbatch](/docs/userguide.html) 
 command.



```
sbatch --mem=10g --cpus-per-task=40 exomiser.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. exomiser.swarm). For example:



```
cd dir1;java -Xms2g -Xmx10g -jar $EXOMISERJAR --analysis examples/test-analysis-exome.yml
cd dir2;java -Xms2g -Xmx10g -jar $EXOMISERJAR --analysis examples/test-analysis-exome.yml
cd dir3;java -Xms2g -Xmx10g -jar $EXOMISERJAR --analysis examples/test-analysis-exome.yml
cd dir4;java -Xms2g -Xmx10g -jar $EXOMISERJAR --analysis examples/test-analysis-exome.yml
cd dir5;java -Xms2g -Xmx10g -jar $EXOMISERJAR --analysis examples/test-analysis-exome.yml
```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f exomiser.swarm -g 10 -t 40 --module exomiser
```

where
 

|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in 
 the swarm command file) 
 | -t *#*  Number of threads required for each process (usually use 40)
 | --module exomiser Loads the vcftools module for each subjob in the swarm 
  | |
 | |
 | |












