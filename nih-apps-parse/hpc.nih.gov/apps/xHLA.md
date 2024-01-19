

document.querySelector('title').textContent = 'Fast and accurate HLA typing from short read sequence data on Biowulf';
**Fast and accurate HLA typing from short read sequence data on Biowulf**


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



The HLA gene complex on human chromosome 6 is one of the
most polymorphic regions in the human genome and contributes
in large part to the diversity of the immune system. Accurate typing
of HLA genes with short-read sequencing data has historically
been difficult due to the sequence similarity between the
polymorphic alleles. xHLA iteratively refines the mapping results at the amino acid level to
achieve high typing accuracy for both class I and II HLA genes.



### References:


* Chao Xiea, Zhen Xuan Yeo, Marie Wong, Jason Piper, Tao Long, Ewen F. Kirkness, William H. Biggs, Ken Bloom,
Stephen Spellman, Cynthia Vierra-Green, Colleen Brady, Richard H. Scheuermann, Amalio Telenti, Sally Howard,
Suzanne Brewerton, Yaron Turpaz, and J. Craig Venter,
Fast and accurate HLA typing from short-read
next-generation sequence data with xHLA, Proc. Natl. Acad. Sci. USA, 2017, vol. 114, N30, p. 8059–8064


Documentation
* [xHLA Github page](https://github.com/humanlongevity/HLA)


Important Notes
* Module Name: xHLA (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Multithreaded
* Unusual environment variables set
	+ **XHLA\_HOME**  xHLA installation directory
	+ **XHLA\_BIN**       xHLA executable directory
	+ **XHLA\_TEST**    xHLA example directory
	+ **XHLA\_DATA**    xHLA data directory* Example file in **$XHLA\_TEST**



**IMPORTANT NOTE:** Before submitting a job running xHLA software to cluster, user(s) should estimate the amount of memory required for running the job by benchmarking their data.
Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=16g --cpus-per-task=14 --gres=lscratch:40**
[user@cn3316 ~]$ **module load xHLA**
[user@cn3316 ~]$ **git clone https://github.com/humanlongevity/HLA.git**
[user@cn3316 ~]$ **mkdir test\_run**
[user@cn3316 ~]$ **xhla --sample\_id test --input\_bam\_path HLA/tests/test.bam --output\_path test\_run**
[07/Feb/2018 15:09:25] INFO - Xie Chao's HLA typing algorithm
[07/Feb/2018 15:09:25] INFO - Sample_id: test Input file: HLA/tests/test.bam
typer.sh parameters: DELETE=false FULL=false
Extracting reads from S3
...
done lpsolve
 [1] "A*01:01"    "A*02:01"    "B*13:02"    "B*37:01"    "C*06:02"   
 [6] "DPB1*04:01" "DQB1*02:01" "DQB1*05:01" "DRB1*07:01" "DRB1*10:01"
pulling non-core exons in
  DQB1*02:01   DQB1*02:01 
"DQB1*02:02" "DQB1*02:12" 
named character(0)
refining solution
        allele rank
 1:    A*01:01  6.0
 2:    A*02:01  3.0
 3:    B*13:02 48.5
 4:    B*37:01 96.0
 5:    C*06:02 10.0
 6: DPB1*04:01   NA
 7: DQB1*02:02   NA
 8: DQB1*05:01 12.0
 9: DRB1*07:01  9.0
10: DRB1*10:01 72.5
     allele rank field1 field2
1:  A*01:01    6      1      1
2: A*01:01L   NA      1      1
3: A*01:01N   NA      1      1
4:  A*01:32   NA      1     32
5:  A*01:45   NA      1     45
6: A*01:56N   NA      1     56
7: A*01:103   NA      1    103
8: A*01:177   NA      1    177
...
Reporting
[07/Feb/2018 15:18:50] INFO - Successfully wrote output file
[07/Feb/2018 15:18:50] INFO - HLA typing: shutting down.
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```

Note that the command **xhla** will accept both the relative and absolute paths to the input data.

Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. xhla.sh). For example:



```

#!/bin/bash
module load xHLA
mkdir HLA/test
xhla --sample_id test --input_bam_path HLA/tests/test.bam --output_path HLA/test

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
**sbatch [--cpus-per-task=#] [--mem=#] xHLA.sh**
```







