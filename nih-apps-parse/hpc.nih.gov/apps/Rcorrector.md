

document.querySelector('title').textContent = 'Rcorrector: Error correction for Illumina RNA-seq reads ';
**Rcorrector: Error correction for Illumina RNA-seq reads** 


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



Rcorrector implements a k-mer based method to correct random sequencing errors 
in Illumina RNA-seq reads. Rcorrector uses a De Bruijn graph to compactly represent 
all trusted k-mers in the input reads. Unlike WGS read correctors, 
which use a global threshold to determine trusted k-mers, Rcorrector computes a local
threshold at every position in a read.



### References:


* Li Song and Liliana Florea   

 *Rcorrector: efficient and accurate error
correction for Illumina RNA-seq reads*  

[GigaScience](https://gigascience.biomedcentral.com/articles/10.1186/s13742-015-0089-y)  2015, **4**, 48


Documentation
* [Rcorrector on GitHub](https://github.com/mourisl/Rcorrector)


Important Notes
* Module Name: Rcorrector (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set
	+ **RCORRECTOR\_HOME**  installation directory
	+ **RCORRECTOR\_BIN**       executable directory
	+ **RCORRECTOR\_DATA**  sample data dorectory
	 + **RCORRECTOR\_SRC**       source code directory
	 + **RCORRECTOR\_DOC**       documentation directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=4g**
[user@@cn3200 ~]$**module load Rcorrector** 
[+] Loading gcc  7.2.0  ... 
[+] Loading jellyfish  2.2.7 
[+] Loading Rcorrector 1.0.3.1  ...
[user@biowulf]$ **rcorrector**
Usage: ./rcorrector [OPTIONS]
OPTIONS:
Required parameters:
	-r seq_file: seq_file is the path to the sequence file. Can use multiple -r to specifiy multiple sequence files
	-p seq_file_left seq_file_right: the paths to the paired-end data set. Can use multiple -p to specifiy multiple sequence files
	-i seq_file: seq_file is the path to the interleaved mate-pair sequence file. Can use multiple -i
	-c jf_dump: the kmer counts dumped by JellyFish
	-k kmer_length
Other parameters:
	-od output_file_directory (default: ./)
	-t number of threads to use (default: 1)
	-maxcor INT: the maximum number of correction every 100bp (default: 8)
	-maxcorK INT: the maximum number of correction within k-bp window (default: 4)
	-wk FLOAT: the proportion of kmers that are used to estimate weak kmer count threshold (default: 0.95)
	-stdout: output the corrected sequences to stdout (default: not used)
	-verbose: output some correction information to stdout (default: not used)
[user@@cn3200 ~]$ **run\_rcorrector.pl** 
Usage: perl ./run_rcorrector.pl [OPTIONS]
OPTIONS:
Required parameters:
	-s seq_files: comma separated files for single-end data sets
	-1 seq_files_left: comma separated files for the first mate in the paried-end data sets
	-2 seq_files_right: comma separated files for the second mate in the paired-end data sets
	-i seq_files_interleaved: comma sperated files for interleaved paired-end data sets
Other parameters:
	-k kmer_length (<=32, default: 23)
	-od output_file_directory (default: ./)
	-t number_of_threads (default: 1)
	-maxcorK INT: the maximum number of correction within k-bp window (default: 4)
	-wk FLOAT: the proportion of kmers that are used to estimate weak kmer count threshold, lower for more divergent genome (default: 0.95)
	-ek expected_number_of_kmers: does not affect the correctness of program but affect the memory usage (default: 100000000)
	-stdout: output the corrected reads to stdout (default: not used)
	-verbose: output some correction information to stdout (default: not used)
	-stage INT: start from which stage (default: 0)
		0-start from begining(storing kmers in bloom filter);
		1-start from count kmers showed up in bloom filter;
		2-start from dumping kmer counts into a jf_dump file;
		3-start from error correction.
[user@@cn3200 ~]$ **run\_rcorrector.pl -1 $RCORRECTOR\_DATA/sample\_read1.fq -2 $RCORRECTOR\_DATA/sample\_read2.fq**
RPut the kmers into bloom filter
jellyfish bc -m 23 -s 100000000 -C -t 1 -o tmp_ffad73bdce91fd172154077d85bac6cc.bc /usr/local/apps/Rcorrector/1.0.3.1/sample_data/sample_read1.fq /usr/local/apps/Rcorrector/1.0.3.1/sample_data/sample_read2.fq 
Count the kmers in the bloom filter
jellyfish count -m 23 -s 100000 -C -t 1 --bc tmp_ffad73bdce91fd172154077d85bac6cc.bc -o tmp_ffad73bdce91fd172154077d85bac6cc.mer_counts /usr/local/apps/Rcorrector/1.0.3.1/sample_data/sample_read1.fq /usr/local/apps/Rcorrector/1.0.3.1/sample_data/sample_read2.fq 
Dump the kmers
jellyfish dump -L 2 tmp_ffad73bdce91fd172154077d85bac6cc.mer_counts > tmp_ffad73bdce91fd172154077d85bac6cc.jf_dump
Error correction
/usr/local/apps/Rcorrector/1.0.3.1/bin/rcorrector   -p /usr/local/apps/Rcorrector/1.0.3.1/sample_data/sample_read1.fq /usr/local/apps/Rcorrector/1.0.3.1/sample_data/sample_read2.fq -c tmp_ffad73bdce91fd172154077d85bac6cc.jf_dump
Stored 253 kmers
Weak kmer threshold rate: 0.010000 (estimated from 0.950/1 of the chosen kmers)
Bad quality threshold is '!'
Processed 70 reads
	Corrected 21 bases.


```

End the interactive session:

```

[user@cn3200 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```

Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. rcorrector.sh). For example:



```

#!/bin/bash
#SBATCH --mem=4g
module load Rcorrector           
run_rcorrector.pl -1 $RCORRECTOR_DATA/sample_read1.fq -2 $RCORRECTOR_DATA/sample_read2.fq

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch rcorrector.sh 
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. rcorrector.swarm). For example:



```

#!/bin/bash
run_rcorrector.pl -1 $RCORRECTOR_DATA/sample_read1.fq -2 $RCORRECTOR_DATA/sample_read2.fq

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f rcorrector.swarm  -g 4 
```





