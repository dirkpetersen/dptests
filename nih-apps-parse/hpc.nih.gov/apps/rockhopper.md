

document.querySelector('title').textContent = 'Rockhopper on Biowulf';
Rockhopper on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
 |




From the Rockhopper manual:




>  Rockhopper is a comprehensive and user-friendly system for
>  computational analysis of bacterial RNA-seq data. As input, Rockhopper
>  takes RNA sequencing reads output by high-throughput sequencing technology
>  (FASTQ, QSEQ, FASTA, SAM, or BAM files).
> 


Rockhopper has a graphical interface which is described in the Rockhopper manual.
For the cluster, the most relevant mode of usage is the command line interface.


### References:


* B. Tjaden. *De novo assembly of bacterial transcriptomes from RNA-seq data.*
Genome Biology 16:1 (2015)
[PubMed](https://www.ncbi.nlm.nih.gov/pubmed/25583448) | 
[PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4316799/) | 
[Journal](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0572-2)



Documentation
* Rockhopper[home page](https://cs.wellesley.edu/~btjaden/Rockhopper/index.html)


Important Notes
* Module Name: rockhopper (see [the modules page](/apps/modules.html) for more information)
* Rockhopper is a multithreaded application (`-p`). Please match the number of allocated 
 CPUs with the number of threads
* The `$ROCKHOPPER_JAR` variable is set to the path of the jar file
* Example files in `$ROCKHOPPER_TEST_DATA`



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and
run the program. If using the GUI, make sure that you have a graphical
connection to biowulf (NX or X11 forwarding though ssh). Sample session:



```

[user@biowulf]$ **sinteractive --mem=5g --cpus-per-task=4 --gres=lscratch:20**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144]$ **module load rockhopper**

[user@cn3144]$ # start up the GUI. Need to know the http proxy host for this.
[user@cn3144]$ # Note that the proxy may change from session to session
[user@cn3144]$ **echo $http\_proxy**
http://dtn03-e0:3128
[user@cn3144]$ **java -Dhttp.proxyHost=dtn03-e0 -Dhttp.proxyPort=3128 -jar $ROCKHOPPER\_JAR**

```

For more details on using the GUI, please see the Rockhopper manual. Now let's use
the command line interface to do a reference based analysis of a small *Mycoplasma genitalium*
data set:



```

[user@cn3144]$ **cd /lscratch/$SLURM\_JOB\_ID**
[user@cn3144]$ **cp -r ${ROCKHOPPER\_TEST\_DATA:-none}/\* .**
[user@cn3144]$ **tree**
.
|-- [user    16M]  Example_Condition1.fastq
|-- [user    16M]  Example_Condition2.fastq
`-- [user   4.0K]  Mycoplasma_genitalium_G37
    |-- [user   575K]  NC_000908.fna
    |-- [user     46]  NC_000908.fna.fai
    |-- [user    10K]  NC_000908.genome
    |-- [user    59K]  NC_000908.gff
    |-- [user    37K]  NC_000908.ptt
    `-- [user   2.3K]  NC_000908.rnt

[user@cn3144]$ **mkdir tmp**
[user@cn3144]$ **java -Xmx4000m -Djava.io.tmpdir=$PWD/tmp -cp $ROCKHOPPER\_JAR Rockhopper \
 -g $PWD/Mycoplasma\_genitalium\_G37 \
 -p $SLURM\_CPUS\_PER\_TASK \
 -L cond1,cond2 \
 -o $PWD/results \
 Example\_Condition1.fastq Example\_Condition1.fastq**
Aligning sequencing reads from file:    Example_Condition1.fastq
Total reads:                    137473
Successfully aligned reads:     125676  91%     (Mycoplasma genitalium G37 chromosome)
        Aligning (sense) to protein-coding genes:       96%
        Aligning (antisense) to protein-coding genes:   0%
        Aligning (sense) to ribosomal RNAs:     0%
        Aligning (antisense) to ribosomal RNAs: 0%
        Aligning (sense) to transfer RNAs:      1%
[...snip...]

[user@cn3144]$ **tree results**
results/
|-- [user   4.0K]  genomeBrowserFiles
|   |-- [user   1.4M]  Example_Condition1_NC_000908_v1_m15_aT_d500_l33_fr_cF.minus.wig
|   |-- [user   1.5M]  Example_Condition1_NC_000908_v1_m15_aT_d500_l33_fr_cF.plus.wig
|   |-- [user   1.1M]  NC_000908_diffExpressedGenes.wig
|   |-- [user   1.1M]  NC_000908_ncRNAs.wig
|   |-- [user   1.3M]  NC_000908_operons.wig
|   `-- [user   1.1M]  NC_000908_UTRs.wig
|-- [user   4.0K]  intermediary
|   `-- [user   186K]  Example_Condition1_NC_000908_v1_m15_aT_d500_l33_fr_cF.gz
|-- [user   4.4K]  NC_000908_operons.txt
|-- [user    37K]  NC_000908_transcripts.txt
`-- [user   2.0K]  summary.txt

[user@cn3144]$ **cp -r results /data/$USER**

[user@cn3144]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. rockhopper.sh), which uses the input file 'rockhopper.in'. For example:



```

#!/bin/bash
wd=$PWD

module load rockhopper/2.0.3
cd /lscratch/$SLURM_JOB_ID
mkdir tmp

cp -r ${ROCKHOPPER_TEST_DATA:-none}/* .
java -Xmx4000m -Djava.io.tmpdir=$PWD/tmp -cp $ROCKHOPPER_JAR Rockhopper \
    -g $PWD/Mycoplasma_genitalium_G37 \
    -p $SLURM_CPUS_PER_TASK \
    -L cond1,cond2 \
    -o $PWD/results \
    Example_Condition1.fastq Example_Condition1.fastq
cp -r results $wd

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=4 --mem=5g rockhopper.sh
```







