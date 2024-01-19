

document.querySelector('title').textContent = 'RSEM on Biowulf';
RSEM on Biowulf


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

  RSEM is a software package for estimating gene and isoform expression levels 
 from RNA-Seq data. The RSEM package provides an user-friendly interface, 
 supports threads for parallel computation of the EM algorithm, single-end 
 and paired-end read data, quality scores, variable-length reads and RSPD 
 estimation. In addition, it provides posterior mean and 95% credibility 
 interval estimates for expression levels. For visualization, It can generate 
 BAM and Wiggle files in both transcript-coordinate and genomic-coordinate. 
 Genomic-coordinate files can be visualized by both UCSC Genome browser and 
 Broad Institute’s Integrative Genomics Viewer (IGV). Transcript-coordinate 
 files can be visualized by IGV. RSEM also has its own scripts to generate 
 transcript read depth plots in pdf format. The unique feature of RSEM is, 
 the read depth plots can be stacked, with read depth contributed to unique 
 reads shown in black and contributed to multi-reads shown in red. In addition, 
 models learned from data can also be visualized. Last but not least, RSEM 
 contains a simulator.


### References:

 * [Bo Li](http://bli25ucb.github.io/) and [Colin 
 Dewey](https://www.biostat.wisc.edu/%7Ecdewey/) designed the RSEM algorithm. Bo Li implemented the RSEM software. 
 [Peng Liu](https://www.biostat.wisc.edu/%7Ecdewey/group.html) 
 contributed the STAR aligner options


Documentation
* <https://github.com/bli25ucb/RSEM_tutorial>


Important Notes
December 2020: While RSEM v1.3.3 is available as a module, version 1.3.2 is set as default because of a bug in version 1.3.3, documented [here](https://github.com/deweylab/RSEM/issues/152). Running rsem-calculate-expression with temporary directory set to /lscratch (as recommended below) throws an error. As a result, we suggest using RSEM v1.3.2.


* Module Name: rsem (see [the modules 
 page](/apps/modules.html) for more information) 
 * Multithreaded and singlethreaded
 * Rsem has shown to use quite a bit of temporary directory space. Please 
 make sure to use "-temporary-folder" flag in "rsem-run-prsem-testing-procedure" 
 and "rsem-calculate-expression" when submit batch job. See example below. 
 * Reference data in /fdb/rsem/



Based on the tutorial <https://github.com/bli25ucb/RSEM_tutorial>


The following steps were already finished on our systems. User can use the index files in the shared area. Some of the steps can not be run by users due to permission setting.


Biowulf > $ mkdir /fdb/rsem  ## users do now have permission to do the following steps in /fdb/rsem but can create directory under personal space such as /data/$USER and replace the following /fdb/rsem with the directory created.


 Biowulf > $ sinteractive --mem=10g  

 cn1234 > $ module load rsem  

 cn1234 > $ cd /fdb/rsem  

 cn1234 > $ ln -s /fdb/ensembl/pub/release-82/fasta/mus\_musculus/dna/Mus\_musculus.GRCm38.dna.toplevel.fa  Mus\_musculus.GRCm38.dna.toplevel.fa  

 cn1234 > $ ln -s /fdb/ensembl/pub/release-82/gtf/mus\_musculus/Mus\_musculus.GRCm38.82.chr.gtf  Mus\_musculus.GRCm38.82.chr.gtf


Build reference from genome (Users can use the reference directly without building) :  

 cn1234 > $ rsem-prepare-reference --bowtie2 --gtf Mus\_musculus.GRCm38.82.chr.gtf Mus\_musculus.GRCm38.dna.toplevel.fa ref\_from\_genome/mouse\_ref  



Build reference from ensemble transcripts (Users can use the reference directly without building) :  

 Downloade mouse\_ref\_building\_from\_transcripts.tar.gz from https://www.dropbox.com/s/ie67okalzaw8zzj/mouse\_ref\_building\_from\_transcripts.tar.gz?dl=0  

 cn1234 > $ tar xvfz mouse\_ref\_building\_from\_transcripts.tar.gz  

 This created two files: mouse\_ref.fa and mouse\_ref\_mapping.txt


cn1234 > $ rsem-prepare-reference --transcript-to-gene-map mouse\_ref\_mapping.txt --bowtie2 mouse\_ref.fa ref\_from\_transcripts/mouse\_ref


 

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

[user@cn3144 ~]$ **module load rsem** bowtie STAR
[user@cn3144 ~]$ **rsem-calculate-expression -p 2 --paired-end --bowtie2 \
 --estimate-rspd --append-names --output-genome-bam \
 /data/user/rsem/SRR937564\_1.fastq /data/user/rsem/SRR937564\_2.fastq \
 /fdb/rsem/ref\_from\_genome/mouse\_ref\_CRCh38.82 LPS\_6h**
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. rsem.sh). For example:



```

#!/bin/bash

module load rsem  bowtie  STAR
cd /data/$USER/dir
rsem-calculate-expression -p $SLURM_CPUS_PER_TASK --temporary-folder /lscratch/$SLURM_JOBID/tempdir \
--paired-end --bowtie2 --estimate-rspd --append-names --output-genome-bam SRR937564_1.fastq SRR937564_2.fastq \
/fdb/rsem/ref_from_genome/mouse_ref LPS_6h

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=4 --mem=10g --gres=lscratch:100 rsem.sh
```

For multiple threaded job, such as -p used in rsem-calculate-expression, 
 use '--cpus-per-task' on the command line and the $SLURM\_CPUS\_PER\_TASK will 
 be assigned automatically to the number user assigned (4 in this case).


 To use temporary local space, use "--gres=lscratch:XXX" which will allocate 
 XXX gb of space for the job. The space can be accessed by using "/lscratch/$SLURM\_JOBID/directoryNameYouChoose" 
 as shown in the script. 


For more memory requirement (default 2xcpus=8gb in this case), use --mem 
 flag


 


Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources. Create a swarmfile (e.g. rsem.swarm). For example:



```


  cd /data/$USER/dir1; rsem command1; rsem command2
  cd /data/$USER/dir2; rsem command1; rsem command2
  cd /data/$USER/dir3; rsem command1; rsem command2
```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f swarmfile -t 4 --module rsem,bowtie,STAR
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module rsem Loads the rsem module for each subjob in the swarm  | |
 | |
 | |










