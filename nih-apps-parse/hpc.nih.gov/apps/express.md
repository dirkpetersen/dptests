

document.querySelector('title').textContent = ' Express on Biowulf';

Express on Biowulf



|  |
| --- |
| 
Quick Links
[Running an interactive job on Biowulf](#int)
[Single Express job on Biowulf](#serial)
[Swarm of Express jobs](#swarm)
[Documentation](#doc)
 |



eXpress is a streaming tool for quantifying the abundances of a set of target sequences from sampled subsequences. Example applications include 
transcript-level RNA-Seq quantification, allele-specific/haplotype expression analysis (from RNA-Seq), transcription factor binding quantification in ChIP-Seq, 
and analysis of metagenomic data. It is based on an online-EM algorithm [[more info here](http://bio.math.berkeley.edu/eXpress/overview.html#ref1)] that 
results in space (memory) requirements proportional to the total size of the target sequences and time requirements that are proportional to the number of sampled 
fragments. Thus, in applications such as RNA-Seq, eXpress can accurately quantify much larger samples than other currently available tools greatly reducing 
computing infrastructure requirements. eXpress can be used to build lightweight high-throughput sequencing processing pipelines when coupled with a streaming 
aligner (such as Bowtie), as output can be piped directly into eXpress, effectively eliminating the need to store read alignments in memory or on disk.



The environment variable(s) need to be set properly first. The easiest way to do this is by using the 
 [modules](modules.html) commands, '**module load express**', as in the example below. 


```

biowulf% module load express

```

This command will load the required `bowtie/1` and `samtools` as well as `express`.


Express was developed by Adam Roberts and Lior Pachter at the University of California, Berkeley. [[Express website](http://bio.math.berkeley.edu/eXpress/index.html)]


Running an interactive job on Biowulf

The sample sessions in this document use the example sets from the Express documentation. In the following two sub-sections, you will run eXpress on a sample RNA-Seq dataset 
with simulated reads from UGT3A2 and the HOXC cluster using the human genome build hg18. Both the transcript sequences (transcripts.fasta) and raw reads (reads\_1.fastq, 
reads\_2.fastq) can be found in the `/usr/local/apps/express/sample_data` directory. 

Before you begin, you must prepare your Bowtie index. Since you wish to allow many multi-mappings, it is useful to build the index with a small offrate (in this case 1). 
The smaller the offrate, the larger the index and the faster the mapping. If you have disk space to spare, always use an offrate of 1.

1. Build the index with the following commands.

```

$ **sinteractive --gres=lscratch:40**
$ **cp -rp /usr/local/apps/express/sample\_data /lscratch/$SLURM\_JOB\_ID/**
$ **cd /lscratch/$SLURM\_JOB\_ID/sample\_data**
$ **module load express**
$ **bowtie-build --offrate 1 transcripts.fasta transcript**

```


This command will populate your directory with several index files that allow Bowtie to more easily align reads to the transcripts.

You can now map the reads to the transcript sequences using the following Bowtie command, which outputs in SAM (-S), allows for unlimited multi-mappings (-a), 
a maximum insert distance of 800 bp between the paired-ends (-X 800), and 3 mismatches (-v 3). The first three options (a,S,X) are highly recommended for best results. 
You should also allow for many mismatches, since eXpress models sequencing errors. Furthermore, you will want to take advantage of multiple processors when mapping large files 
using the -poption. See the Bowtie Manual for more details on various parameters and options. 

- The SAM output from Bowtie is piped into SAMtools in order to compress it to BAM format. This conversion is optional, but will greatly reduce the size of the alignment file.


```

$ **bowtie -aS -X 800 --offrate 1 -v 3 transcript -1 reads\_1.fastq -2 reads\_2.fastq | samtools view -Sb - > hits.bam**

```


- Once you have aligned your reads to the transcriptome and stored them in a SAM or BAM file, you can run eXpress in default mode with the command:


```

$ **module load express**
$ **express transcripts.fasta hits.bam**

```


- If you do not wish to store an intermediate SAM/BAM file, you can pipe the Bowtie output directly into eXpress with the command:


```

$ **module load express**
$ **bowtie -aS -X 800 --offrate 1 -v 3 transcript -1 reads\_1.fastq -2 reads\_2.fastq | express transcripts.fasta** 
$ **exit** 
  
```



Running a single Express batch job on Biowulf

The following batch script uses the same commands as in the preceding example.


```

#!/bin/bash
#
# this file is called express.bat

module load express
cd /data/$USER/express
cp -rp $EXPRESS_ROOT/sample_data  .
cd sample_data
bowtie-build --offrate 1 transcripts.fasta transcript
bowtie -p $SLURM_CPUS_PER_TASK -aS -X 800 --offrate 1 -v 3 transcript -1 reads_1.fastq -2 reads_2.fastq \
          | express transcripts.fasta 

```


Submit this job with:


```

sbatch  --cpus-per-task=4  express.bat

```


Bowtie is a multi-threaded program, and so the job is submitted to 4 CPUs ('--cpus-per-task=4' in the sbatch command above). 
In the batch script express.bat, bowtie is run with '-p $SLURM\_CPUS\_PER\_TASK' to utilize all 4 cpus. 


Running a swarm of Express batch jobs on Biowulf

Set up a swarm command file along the following lines:

```

bowtie -p $SLURM_CPUS_PER_TASK -aS -X 800 --offrate 1 -v 3 transcript -1 reads_1.fastq -2 reads_2.fastq  \
        | express transcripts.fasta
bowtie -p $SLURM_CPUS_PER_TASK -aS -X 800 --offrate 1 -v 3 transcript -1 reads_3.fastq -2 reads_4.fastq \ 
       | express transcripts.fasta
etc...

```

Submit this job with:

```

swarm -t 4 -f cmdfile --module express

```


The parameter '-t 4' tells swarm to allocate 4 CPUs for each command above. Express is single-threaded, but bowtie will make
use of the 4 CPUS ('-p $SLURM\_CPUS\_PER\_TASK' in the swarm command file above).

By default, each line in the swarm command file above will be executed on 4 CPU and can utilize up to 8 GB of memory.
If your commands require more than 8 GB, you should submit with:

```

swarm -g # -f cmdfile --module express

```

where '#' is the number of GigaBytes of memory required for a single process (1 line in your swarm command file).

Documentation

<http://bio.math.berkeley.edu/eXpress/manual.html>
















































