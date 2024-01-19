

document.querySelector('title').textContent = ' Transvar on Biowulf';
 Transvar on Biowulf


|  |
| --- |
| 
Quick Links
[Interactive Job on Biowulf](#int)
[Single Batch Job on Biowulf](#serial)
[Swarm of Jobs](#swarm)
[Documentation](#doc)
 |


**TransVar** is a versatile annotator for 3-way conversion and annotation among genomic characterization(s) of mutations (e.g., chr3:g.178936091G>A) and transcript-dependent annotation(s) (e.g., PIK3CA:p.E545K or PIK3CA:c.1633G>A, or NM\_006218.2:p.E545K, or NP\_006266.2:p.G240Afs\*50). It is particularly designed with the functionality of resolving ambiguous mutation annotations arising from differential transcript usage. TransVar keeps awareness of the underlying unknown transcript structure (exon boundary, reference amino acid/base) while performing reverse annotation (via fuzzy matching from protein level to cDNA level). TransVar has the following features:


* supports HGVS nomenclature
* supports input from gene name, transcript ID, protein ID, UniProt ID and other aliases.
* supports both left-alignment and right-alignment convention in reporting indels and duplications.
* supports annotation of a region based on a transcript-dependent characterization
* supports mutations at both coding region and intronic/UTR regions
* supports noncoding RNA annotation
* supports VCF inputs
* supports long haplotype decomposition
* supports single nucleotide variation (SNV), insertions and deletions (indels) and block substitutions
* supports transcript annotation from commonly-used databases such as Ensembl, NCBI RefSeq and GENCODE etc
* supports GRCh36, 37, 38 (human), GRCm38 (mouse), NCBIM37 (mouse)
* supports >60 other genomes available from Ensembl
* functionality of forward annotation.


**Note 1:**


Transvar Hg19 files are located under /fdb/transvar
**Note 2:** 
To avoid /home/$USER disk quota filled up, create a link to point '/home/$USER/.transvar.download' to '/data/$USER/transvar' first.
**Note 3:** 
In version 2.5.9 and newer, download of new/current versions of the annotations is broken. Older versions can be downloaded either manually or by temporarily using an older version, such as 2.4, to access them from the old web site.


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

[user@cn3144 ~]$ **module load transvar**

[user@cn3144 ~]$ **transvar config --download\_anno --refversion hg19**

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Submitting a single batch job
1. Create a script file. The file will contain the lines similar to the
 lines below. Modify the path of program location before running.




```
#!/bin/bash 

module load transvar
cd /data/$USER/somewhere
transvar config --download_anno --refversion hg19
....
....
```

2. Submit the script on Biowulf. 



$ sbatch myscript





Submitting a swarm of jobs
Using the 'swarm' utility, one can submit many jobs to the cluster to run concurrently. 


Set up a swarm command file (eg /data/$USER/cmdfile). Here is a sample file:



```
cd /data/user/run1/; transvar config --download_anno --refversion hg19
cd /data/user/run2/; transvar config --download_anno --refversion hg19
cd /data/user/run3/; transvar config --download_anno --refversion hg19
........

```

The **-f** flag is required to specify swarm file name.  

Submit the swarm job:

```
$ swarm -f swarmfile --module transvar
```

- Use -g flag for more memory requirement (default 1.5gb per line in swarmfile)


For more information regarding running swarm, see [swarm.html](http://hpc.cit.nih.gov/apps/swarm.html)


 



Documentation
<https://bitbucket.org/wanding/transvar>
























