

document.querySelector('title').textContent = 'Gemini on Biowulf';
Gemini on Biowulf


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



GEMINI (GEnome MINIng) is designed to be a flexible framework for exploring genetic variation in the context of the wealth of genome annotations available for the human genome. By placing genetic variants, sample genotypes, and useful genome annotations into an integrated database framework, GEMINI provides a simple, flexible, yet very powerful system for exploring genetic variation for for disease and population genetics. 


### References:


* Paila U, Chapman BA, Kirchner R, Quinlan AR (2013) GEMINI: Integrative Exploration of Genetic Variation and Genome Annotations. PLoS Comput Biol 9(7): e1003153. doi:10.1371/journal.pcbi.1003153


Documentation
* [Gemini Main Site](https://gemini.readthedocs.io/en/latest/)


Important Notes
* Module Name: gemini (see [the modules page](/apps/modules.html) for more information)
* Singlethreaded
* Example files in /fdb/gemini* Reference data in /fdb/gemini/



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

[user@cn3144 ~]$ **cd /data/$USER/gemini**

[user@cn3144 ~]$ **cp -r /fdb/gemini/gemini/test .**

[user@cn3144 ~]$ **cp /fdb/gemini/gemini/master-test.sh .**

[user@cn3144 ~]$ **module load gemini**

[user@cn3144 ~]$ **bash master-test.sh**
Bgzipping test.query.vcf into test.query.vcf.gz.
Indexing test.query.vcf.gz with grabix.
Loading 879 variants.
Breaking test.query.vcf.gz into 2 chunks.
Loading chunk 0.
Loading chunk 1.
Done loading 879 variants in 1 chunks.
    region.t01...\c
ok
    region.t02...\c
ok
    region.t03...\c
ok
    region.t04...\c
[...]

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. gemini.sh). For example:



```

#!/bin/bash
cd /data/$USER/gemini
cp -r /fdb/gemini/gemini/test .
cp /fdb/gemini/gemini/master-test.sh .
bash master-test.sh

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --mem=5g gemini.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.






