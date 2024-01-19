

document.querySelector('title').textContent = 'locuszoom on Biowulf';
locuszoom on Biowulf


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



LocusZoom is a tool for visualizing results of genome wide association studies
at an individual locus along with other relevant information like gene models,
linkage disequillibrium coefficients, and estimated local recombination rates.



LocusZoom uses association results in METAL or EPACTS formatted files along
with it's own source of supporting data (see below) to generate graphs.


### References:


* Randall J. Pruim et al.*LocusZoom: regional visualization of 
 genome-wide association scan results*. Bioinformatics 2010, 26:2336-2337.
 [Pubmed](http://www.ncbi.nlm.nih.gov/pubmed/20634204) | 
 [PMC](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2935401/) | 
 [Journal](http://bioinformatics.oxfordjournals.org/content/26/18/2336.abstract?keytype=ref&ijkey=IZbX3XuU1KUZXYZ)


Documentation
* [Home page](http://locuszoom.sph.umich.edu/locuszoom/)
* [Web service](http://locuszoom.sph.umich.edu/locuszoom/genform.php?type=yourdata) [hosted at UMich]
* [Documentation](http://genome.sph.umich.edu/wiki/LocusZoom_Documentation)
* [LocusZoom
 command line documentation](http://genome.sph.umich.edu/wiki/LocusZoom_Standalone#Download)


Important Notes
* Module Name: locuszoom (see [the modules page](/apps/modules.html) for more information)
* Example files in `$LOCUSZOOM_TEST_DATA`
### Supporting data


LocusZoom contains a number of data files used in generating graph annotations:

+ Linkage disequillibrium
+ SNP, gene, exon positions
+ recombination rates

All data for LocusZoom is stored in **/fdb/locuszoom/[version]**. 
LocusZoom has been configured to automatically find all required information.




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

[user@cn3144]$ **module load locuszoom**
[user@cn3144]$ **locuszoom -h**
locuszoom -h
+---------------------------------------------+
| LocusZoom 1.3 (06/20/2014)                  |
| Plot regional association results           |
| from GWA scans or candidate gene studies    |
+---------------------------------------------+

usage: locuszoom [options]

  -h, --help
    show this help message and exit

  --metal <string>
    Metal file.
[...snip...]

```

Draw a diagram of the associations between SNPs and HDL observed in [Kathiresan et al, 2009](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2881676/)
around the FADS1 gene



```

[user@cn3144]$ **locuszoom \
 --metal /usr/local/apps/locuszoom/TEST\_DATA/examples/Kathiresan\_2009\_HDL.txt \
 --refgene FADS1**

```


![LocusZoom output](/images/locuszoom_chr11_61303672-61361105.png)

Exit the interactive session



```

[user@cn3144]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. locuszoom.sh), which uses the input file 'locuszoom.in'. For example:



```

#! /bin/bash
#SBATCH --mail-type=END
# this is locuszoom.sh
set -e
function fail {
  echo "$@" >&2
  exit 1
}

module load locuszoom/1.3 || fail "Could not load locuszoom module"


mf=$LOCUSZOOM_TEST_DATA/examples/Kathiresan_2009_HDL.txt
locuszoom --metal=$mf --refgene FADS1 --pop EUR --build hg19 \
  --source 1000G_March2012 \
  --gwas-cat whole-cat_significant-only

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] locuszoom.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. locuszoom.swarm). For example:



```

locuszoom --metal=$LOCUSZOOM_TEST_DATA/examples/Kathiresan_2009_HDL.txt \
  --refgene FADS1
locuszoom --metal=$LOCUSZOOM_TEST_DATA/examples/Kathiresan_2009_HDL.txt \
  --refgene PLTP
locuszoom --metal=$LOCUSZOOM_TEST_DATA/examples/Kathiresan_2009_HDL.txt \
  --refgene ANGPTL4

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f locuszoom.swarm [-g #] [-t #] --module locuszoom
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g #  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t #  Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module locuszoom  Loads the locuszoom module for each subjob in the swarm 
 | |
 | |
 | |








