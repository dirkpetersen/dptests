

document.querySelector('title').textContent = 'shmlast: an improved implementation of Conditional Reciprocal Best Hits with LAST and Python ';
**shmlast: an improved implementation of Conditional Reciprocal Best Hits with LAST and Python** 


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



shmlast is a reimplementation of the Conditional Reciprocal Best Hits algorithm for finding potential orthologs between a transcriptome and a species-specific protein database. It uses the LAST aligner and the pydata stack to achieve much better performance while staying in the Python ecosystem.



### References:


* Bo Li, Nathanael Fillmore, Yongsheng Bai, Mike Collins, James A. Thomson, Ron Stewart, and Colin N. Dewey.   

 *Evaluation of de novo transcriptome assemblies from RNA-Seq data*  

[Genome Biology](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0553-5)  2014, **15**:553


Documentation
* [shmlast project home page](https://pypi.org/project/shmlast/)
* [shmlast github page](https://github.com/camillescott/shmlast)


Important Notes
* Module Name: shmlast (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set
	+ **SHMLAST\_HOME**  installation directory
	+ **SHMLAST\_BIN**       executable directory
	+ **SHMLAST\_DATA**  sample data directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=4g**
[user@@cn3200 ~]$**module load shmlast** 
[+] Loading singularity  3.10.0  on cn1113
[+] Loading shmlast  1.6
[user@biowulf]$ **shmlast**
usage: shmlast [-h] [--version] {rbl,crbl} ...

shmlast is a reimplementation of the Conditional Reciprocal Best
Hits algorithm for finding potential orthologs between
a transcriptome and a species-specific protein database. It uses the LAST
aligner and the pydata stack to achieve much better performance while staying in the Python ecosystem. 

positional arguments:
  {rbl,crbl}

optional arguments:
  -h, --help  show this help message and exit
  --version   show program's version number and exit
[user@biowulf]$ 
shmlast 1.2.1 -- Camille Scott, 2016
------------------------------------
subcommand: Conditional Reciprocal Best LAST
doit action: run

--- Begin Task Execution ---
. rename:/usr/local/apps/shmlast/1.2.1/sample\_data/test-transcript.fa:
 \* Python: rename\_input
. rename:/usr/local/apps/shmlast/1.2.1/sample\_data/test-protein.fa:
 \* Python: rename\_input
. translate:.test-transcript.fa:
 \* Python: function translate\_fastx
. lastdb:.test-transcript.fa.pep:
 \* Cmd: `/usr/local/apps/shmlast/1.2.1/bin/lastdb -p -w3 .test-transcript.fa.pep .test-transcript.fa.pep`
. lastdb:.test-protein.fa:
 \* Cmd: `/usr/local/apps/shmlast/1.2.1/bin/lastdb -p -w3 .test-protein.fa .test-protein.fa`
. lastal:.test-protein.fa.x.test-transcript.fa.pep.maf:
 \* Cmd: `cat .test-protein.fa | /usr/local/apps/parallel/20171222/bin/parallel --round-robin --pipe -L 2 -N 10000 --gnu -j 1 -a .test-protein.fa /usr/local/apps/shmlast/1.2.1/bin/lastal -D100000.0 .test-transcript.fa.pep > .test-protein.fa.x.test-transcript.fa.pep.maf`
. lastal:.test-transcript.fa.pep.x.test-protein.fa.maf:
 \* Cmd: `cat .test-transcript.fa.pep | /usr/local/apps/parallel/20171222/bin/parallel --round-robin --pipe -L 2 -N 10000 --gnu -j 1 -a .test-transcript.fa.pep /usr/local/apps/shmlast/1.2.1/bin/lastal -D100000.0 .test-protein.fa > .test-transcript.fa.pep.x.test-protein.fa.maf`
. fit\_and\_filter\_crbl\_hits:
 \* Python: do\_crbl\_fit\_and\_filter

```

End the interactive session:

```

[user@cn3200 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```

Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. shmlast.sh). For example:



```

#!/bin/bash
#SBATCH --mem=4g
module load shmlast
cd /data/$USER         
shmlast crbl -q $SHMLAST_DATA/test-transcript.fa -d  $SHMLAST_DATA/test-protein.fa     

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch shmlast.sh 
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. shmlast.swarm). For example:



```

#!/bin/bash
module load shmlast
cd /data/$USER
shmlast crbl -q $SHMLAST_DATA/test-transcript.fa -d  $SHMLAST_DATA/test-protein.fa

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f shmlast.swarm -g 4
```





