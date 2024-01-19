

document.querySelector('title').textContent = 'MAJIQ and VOILA';
MAJIQ and VOILA



|  |
| --- |
| 
Quick Links
[Interactive job](#helix)
[Batch job on Biowulf](#sbatch)
[Swarm of jobs on Biowulf](#swarm)
[Documentation](#doc)
 |



MAJIQ and Voila are two software packages that together detect, quantify, and visualize local splicing variations (LSV) from RNA-Seq data. Conceptually, MAJIQ/Voila can be divided into three modules:
* **MAJIQ Builder**: Uses RNA-Seq (BAM files) and a transcriptome annotation file (GFF3) to define splice graphs and known/novel Local Splice Variations (LSV).
* **MAJIQ Quantifier**: Quantifies relative abundance (PSI) of LSVs and changes in relative LSV abundance (delta PSI) between conditions w/wo replicates.
* **Voila**: A visualization package that combines the output of MAJIQ Builder and MAJIQ Quantifier using interactive D3 components and HTML5. Voila creates interactive summary files with gene splice graphs, LSVs, and their quantification.


**Web sites**
* [Home](http://majiq.biociphers.org/)
* [Biociphers](http://www.biociphers.org/)
* [Discussion Group](https://groups.google.com/forum/#!forum/majiq_voila)


**References**
* [Vaquero-Garcia, Jorge, et al. "A new view of transcriptome complexity and regulation through the lens of local splicing variations." *Elife* 5 (2016): e11752.](https://elifesciences.org/content/5/e11752)
* [Sotillo, Elena, et al. "Convergence of acquired mutations and alternative splicing of CD19 enables resistance to CART-19 immunotherapy." *Cancer discovery* 5.12 (2015): 1282-1295.](http://cancerdiscovery.aacrjournals.org/content/5/12/1282.short)



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
[back to top](majiq.html)  


MAJIQ build jobs require a configuration file. A typical configuration file may look something like this: 

```

[info]
readlen=50
samdir=/data/$USER/my_project
genome=hg19
genome_path=/fdb/igenomes/Homo_sapiens/Ensembl/GRCh37/Sequence/Bowtie2Index/
[experiments]
WT=S1,S2
KO=S4,S5

```


See the [MAJIQ Tutorial](http://majiq.biociphers.org/docs/TUTORIAL%20-%20V1.1.pdf) for more info. 
 
You also need an annotation file in GFF3 format. See more on the [MAJIQ home page](http://majiq.biociphers.org/tools.php).

Once you've written the configuration file and obtained the annotation file, you can start a MAJIQ build session like so:

Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive -c 8 --mem 10g**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load majiq**

[user@cn3144 ~]$ **majiq build ensembl.hg19.gff3 -conf config.file --nthreads 8 --output ./build/**

```


Once MAJIQ build has finished running, MAJIQ psi can be called on the output like so:

```

:[user@cn3144 ~]$ **majiq psi ./build/S1.majiq ./build/S2.majiq --nthreads 8 --output ./psi --name WT**
[user@cn3144 ~]$ **majiq psi ./build/S3.majiq ./build/S4.majiq --nthreads 8 --output ./psi --name KO**

```


Then the deltapsi can be calculated like so:

```

[user@cn3144 ~]$ **majiq deltapsi -grp1 ./build/S1.majiq ./build/S2.majiq --nthreads 8 --output ./dpsi --name WT**
[user@cn3144 ~]$ **majiq deltapsi -grp1 ./build/S3.majiq ./build/S4.majiq --nthreads 8 --output ./dpsi --name KO**

```


When MAJIQ deltapsi has finished, VOILA can be called in turn to visualize the analyzed data. For this step, an [X11 connection](https://hpc.nih.gov/docs/connect.html) will be necessary. See the [MAJIQ Tutorial](http://majiq.biociphers.org/docs/TUTORIAL%20-%20V1.1.pdf) for more details.

Running a single MAJIQ job on Biowulf
[back to top](majiq.html)  


After setting the correct files as described [above](#helix), set up a batch script along the following lines:

```

#!/bin/bash
# this script is called myjob.bat

module load majiq
cd /data/$USER/my_project

majiq build ensembl.hg19.gff3 -conf config.file --nthreads 8 --output ./build/

```



Submit this job with:

```

[user@biowulf ~]$ **sbatch --cpus-per-task=8 myjob.bat**

```


Note that the --cpus-per-task argument in sbatch must match the --nthreads argument in the MAJIQ batch script. Depending on the data, memory requirements may also need to be increased using the --mem or --mem-per-cpu option. See the documentation on [Job Submission](https://hpc.nih.gov/docs/userguide.html#submit) for more info on how to use sbatch.

Running a swarm of MAJIQ jobs on Biowulf
[back to top](majiq.html)  


Set up a swarm command file containing one line for each of your MAJIQ runs. You will also have to generate a separate configuration file for each swarm command so that different experiments are run in different jobs. 
Sample swarm command file



```

majiq build ensembl.hg19.gff3 -conf config.file0 --nthreads 8 --output ./build/
majiq build ensembl.hg19.gff3 -conf config.file1 --nthreads 8 --output ./build/
majiq build ensembl.hg19.gff3 -conf config.file2 --nthreads 8 --output ./build/
majiq build ensembl.hg19.gff3 -conf config.file3 --nthreads 8 --output ./build/

```

Submit this set of runs to the batch system by typing:



```

[user@biowulf ~]$ **swarm --threads-per-process 8 --module majiq -f myjobs.swarm**

```


For details on using swarm see [Swarm on 
Biowulf](/apps/swarm.html).

Documentation
[back to top](majiq.html)  


[Online Tutorial](http://majiq.biociphers.org/docs/TUTORIAL%20-%20V1.1.pdf).




































