

document.querySelector('title').textContent = 'VarScan: Variant calling and somatic mutation/CNV detection for next-generation sequencing data';
**VarScan: Variant calling and somatic mutation/CNV detection for next-generation sequencing data**


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



VarScan is an open source tool for variant detection that is compatible with several short
read aligners. It is capable of detecting SNPs
and indels with high sensitivity and specificity, in both Roche/454
sequencing of individuals and deep Illumina/Solexa sequencing of
pooled samples. VarScan2 detects somatic mutations and copy number alterations (CNAs) in exome data from tumor–normal pairs.



### References:


* Daniel C. Koboldt, Ken Chen, Todd Wylie, David E. Larson, Michael D. McLellan,
Elaine R. Mardis, George M. Weinstock, Richard K. Wilson and Li Ding   

*VarScan: variant detection in massively parallel sequencing of
individual and pooled samples*  
 Bioinformatics, 2009, **25**(17), pp.2283-2285.
* Daniel C. Koboldt, Qunyuan Zhang, David E. Larson, Dong Shen,
Michael D. McLellan, Ling Lin, Christopher A. Miller, Elaine R. Mardis, Li Ding,
and Richard K. Wilson  

*VarScan 2: Somatic mutation and copy number
alteration discovery in cancer by exome sequencing*    
 Genome Res. 2012, **22**, pp.568-576.


Documentation
* [VarScan Github page](https://github.com/dkoboldt/varscan)


Important Notes
* Module Name: VarScan (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set
	+ **VARSCAN\_BIN**       VarScan executable directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive**
[user@cn3316 ~]$ **module load VarScan**

```

Download/prepare sample input data: 

```

[user@cn3316 ~]$ **URL="http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq"**
[user@cn3316 ~]$ **wget $URL/wgEncodeUwRepliSeqBg02esG1bAlnRep1.bam -O myData1.bam**
[user@cn3316 ~]$ **wget $URL/wgEncodeUwRepliSeqBg02esG2AlnRep1.bam -O myData2.bam**
[user@cn3316 ~]$ **samtools sort myData1.bam > myData1\_sorted.bam** 
[user@cn3316 ~]$ **samtools sort myData2.bam > myData2\_sorted.bam** 
[user@cn3316 ~]$ **ln -s /fdb/igenomes/Homo\_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa hg19.fa**
[user@cn3316 ~]$ **samtools mpileup -B -f hg19.fa myData1\_sorted.bam > myData1.pileup**
[user@cn3316 ~]$ **samtools mpileup -B -f hg19.fa myData2\_sorted.bam > myData2.pileup**

```

Run a VarScan executable on the inputs:

```

[user@cn3316 ~]$**varscan somatic myData1.pileup myData2.pileup --output-snp snp --output-indel indel**
Normal Pileup: myData1.pileup
Tumor Pileup: myData2.pileup
NOTICE: While dual input files are still supported, using a single mpileup file (normal-tumor) with the --mpileup 1 setting is strongly recommended.
Min coverage:   8x for Normal, 6x for Tumor
Min reads2:     2
Min strands2:   1
Min var freq:   0.2
Min freq for hom:       0.75
Normal purity:  1.0
Tumor purity:   1.0
Min avg qual:   15
P-value thresh: 0.99
Somatic p-value:        0.05
52831015 positions in tumor
52831015 positions shared in normal
225836 had sufficient coverage for comparison
0 were called Reference
0 were mixed SNP-indel calls and filtered
225836 were called Germline
0 were called LOH
0 were called Somatic
0 were called Unknown
0 were called Variant
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. VarScan.sh). For example:



```

#!/bin/bash
module load VarScan
varscan somatic myData1.pileup myData2.pileup --output-snp snp12 --output-indel indel12
varscan somatic myData3.pileup myData4.pileup --output-snp snp34 --output-indel indel34
...

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] VarScan.sh
```







