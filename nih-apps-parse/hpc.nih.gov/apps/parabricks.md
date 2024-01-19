

document.querySelector('title').textContent = 'parabricks: a software suite for performing secondary analysis of NGS DNA and RNA data ';
**parabricks: a software suite for performing secondary analysis of NGS DNA and RNA data** 


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



The Clara Parabricks toolkit is a set of GPU-accelerated genome analysis tools for secondary 
analysis of next generation sequencing data. It includes GPU-accelerated tools for preprocessing, QC,
alignment, variant calling, and GVCF processing.



### References:


* Mohammed Alser, Zülal Bingöl, et al.  [*Accelerating Genome Analysis: A Primer on an Ongoing Journey.*](https://ieeexplore.ieee.org/abstract/document/9154510)  IEEE Micro, Volume: 40, Issue: 5, 2020.
* Karl R Franke, Erin L Crowgey. [*Accelerating next generation sequencing data analysis: an evaluation of optimized best practices for Genome Analysis Toolkit algorithms*](https://pubmed.ncbi.nlm.nih.gov/32224843/). Genomics & Informatics, Volume 18(1), 2020.


Documentation
* [Parabricks 4.0.0 documentation](https://docs.nvidia.com/clara/parabricks/4.0.0/index.html)


Important Notes
* Module Name: parabricks (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Parabricks requires a GPU with at least CUDA architecture 60. Currently, on the cluster this means
 a V100, V100x, or A100 gpu. You can select from these by specifying --gres=gpu:1 and 
 --constraint="gpua100|gpuv100|gpuv100x".
* Monitor your memory utilization with the [HPC Dashboard](https://hpcnihapps.cit.nih.gov/auth/dashboard/). Running low on memory may cause unpredictable slowdowns in Parabricks software.
* Running out of GPU memory may cause an error like: cudaSafeCall() failed at ParaBricks/src/samGenerator.cu/772: out of memory
+ Some Parabricks programs include a --low-memory option that may help with this

* All tools in Parabricks are used the the pbrun command line tool.
* Unusual environment variables set
	+ **PARABRICKS\_TEST\_DATA**: Parabricks sample data



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
  


```

[user@biowulf]$ **sinteractive --gres=lscratch:40,gpu:1 \
 --mem=32g \
 -c16 \
 --constraint="gpua100|gpuv100|gpuv100x"**

[user@cn4464 ~]$ **module load parabricks** 
[+] Loading parabricks  4.0.0  on cn4464
[+] Loading singularity  3.8.5-1  on cn4464

[user@cn4464 ~]$ **cd /lscratch/$SLURM\_JOB\_ID**

[user@cn4464 ~]$ **tar xf $PARABRICKS\_TEST\_DATA**

[user@cn4464 ~]$ **pbrun --help**
Please visit https://docs.nvidia.com/clara/#parabricks for detailed documentation

usage: pbrun  []
Help: pbrun -h

command can be a TOOL or FULL PIPELINE. Example:
pbrun fq2bam --ref genome.fa --in-fq sample\_1.fq.gz sample\_2.fq.gz --out-bam sample.bam
pbrun germline --ref genome.fa --in-fq sample\_1.fq.gz sample\_2.fq.gz --out-bam sample.bam --out-variants sample.vcf

command options for standalone TOOL
applybqsr - Apply BQSR report to a BAM file and generate a new BAM file
bam2fq - Convert a BAM file to FASTQ
bammetrics - Collect WGS Metrics on a BAM file
bamsort - Sort a BAM file
bqsr - Collect BQSR report on a BAM file
collectmultiplemetrics - Collect multiple classes of metrics on a BAM file
dbsnp - Annotate variants based on a dbsnp
deepvariant - Run GPU-DeepVariant for calling germline variants
fq2bam - Run bwa mem, co-ordinate sorting, marking duplicates, and Base Quality Score Recalibration
genotypegvcf - Convert a GVCF to VCF
haplotypecaller - Run GPU-HaplotypeCaller for calling germline variants
indexgvcf - Index a GVCF file
mutectcaller - Run GPU-Mutect2 for tumor-normal analysis
postpon - Generate the final VCF output of doing mutect pon
prepon - Build an index for PON file, which is the prerequisite to performing mutect pon
rna\_fq2bam - Run RNA-seq data through the fq2bam pipeline
starfusion - Identify candidate fusion transcripts supported by Illumina reads

command options for commonly used FULL PIPELINES
germline - Run the germline pipeline from FASTQ to VCF
deepvariant\_germline - Run the germline pipeline from FASTQ to VCF using a deep neural network analysis
somatic - Run the somatic pipeline from FASTQ to VCF

Information about the software
version - Current version of Parabricks

Please visit https://docs.nvidia.com/clara/#parabricks for detailed documentation

positional arguments:
 command The pipeline or tool to run.

optional arguments:
 -h, --help show this help message and exit

# It may be necessary to add --low-memory on V100 nodes
[user@cn4464 ~]$ **time pbrun fq2bam \
 --ref ./parabricks\_sample/Ref/Homo\_sapiens\_assembly38.fasta \
 --in-fq ./parabricks\_sample/Data/sample\_1.fq.gz ./parabricks\_sample/Data/sample\_2.fq.gz \
 --out-bam fq2bam\_output.bam \
 --tmp-dir /lscratch/$SLURM\_JOB\_ID/tmp**
[Parabricks Options Mesg]: Checking argument compatibility
[Parabricks Options Mesg]: Automatically generating ID prefix
[Parabricks Options Mesg]: Read group created for /lscratch/48192446/parabricks\_sample/Data/sample\_1.fq.gz and
/lscratch/48192446/parabricks\_sample/Data/sample\_2.fq.gz
[Parabricks Options Mesg]: @RG\tID:HK3TJBCX2.1\tLB:lib1\tPL:bar\tSM:sample\tPU:HK3TJBCX2.1
[PB Info 2022-Sep-21 16:11:10] ------------------------------------------------------------------------------
[PB Info 2022-Sep-21 16:11:10] || Parabricks accelerated Genomics Pipeline ||
[PB Info 2022-Sep-21 16:11:10] || Version 4.0.0-1 ||
[PB Info 2022-Sep-21 16:11:10] || GPU-BWA mem, Sorting Phase-I ||
[PB Info 2022-Sep-21 16:11:10] ------------------------------------------------------------------------------
[M::bwa\_idx\_load\_from\_disk] read 0 ALT contigs
[PB Info 2022-Sep-21 16:11:31] GPU-BWA mem
[PB Info 2022-Sep-21 16:11:31] ProgressMeter Reads Base Pairs Aligned
[PB Info 2022-Sep-21 16:12:10] 5043564 580000000
[PB Info 2022-Sep-21 16:12:44] 10087128 1160000000
[PB Info 2022-Sep-21 16:13:16] 15130692 1740000000
[PB Info 2022-Sep-21 16:13:52] 20174256 2320000000
[PB Info 2022-Sep-21 16:14:24] 25217820 2900000000
[PB Info 2022-Sep-21 16:15:00] 30261384 3480000000
[PB Info 2022-Sep-21 16:15:35] 35304948 4060000000
[PB Info 2022-Sep-21 16:16:10] 40348512 4640000000
[PB Info 2022-Sep-21 16:16:43] 45392076 5220000000
[PB Info 2022-Sep-21 16:17:18] 50435640 5800000000
[PB Info 2022-Sep-21 16:17:39]
GPU-BWA Mem time: 368.713119 seconds
[PB Info 2022-Sep-21 16:17:39] GPU-BWA Mem is finished.


[main] CMD: /usr/local/parabricks/binaries//bin/bwa mem -Z ./pbOpts.txt /lscratch/48192446/parabricks\_sample/Ref/Homo\_sapiens\_assembly38.fasta /lscratch/48192446/parabricks\_sample/Data/sample\_1.fq.gz /lscratch/48192446/parabricks\_sample/Data/sample\_2.fq.gz @RG\tID:HK3TJBCX2.1\tLB:lib1\tPL:bar\tSM:sample\tPU:HK3TJBCX2.1
[main] Real time: 389.507 sec; CPU: 5828.472 sec
[PB Info 2022-Sep-21 16:17:39] ------------------------------------------------------------------------------
[PB Info 2022-Sep-21 16:17:39] || Program: GPU-BWA mem, Sorting Phase-I ||
[PB Info 2022-Sep-21 16:17:39] || Version: 4.0.0-1 ||
[PB Info 2022-Sep-21 16:17:39] || Start Time: Wed Sep 21 16:11:10 2022 ||
[PB Info 2022-Sep-21 16:17:39] || End Time: Wed Sep 21 16:17:39 2022 ||
[PB Info 2022-Sep-21 16:17:39] || Total Time: 6 minutes 29 seconds ||
[PB Info 2022-Sep-21 16:17:39] ------------------------------------------------------------------------------
[PB Info 2022-Sep-21 16:17:42] ------------------------------------------------------------------------------
[PB Info 2022-Sep-21 16:17:42] || Parabricks accelerated Genomics Pipeline ||
[PB Info 2022-Sep-21 16:17:42] || Version 4.0.0-1 ||
[PB Info 2022-Sep-21 16:17:42] || Sorting Phase-II ||
[PB Info 2022-Sep-21 16:17:42] ------------------------------------------------------------------------------
[PB Info 2022-Sep-21 16:17:42] progressMeter - Percentage
[PB Info 2022-Sep-21 16:17:42] 0.0 0.00 GB
[PB Info 2022-Sep-21 16:17:52] 81.5 0.00 GB
[PB Info 2022-Sep-21 16:18:02] Sorting and Marking: 20.003 seconds
[PB Info 2022-Sep-21 16:18:02] ------------------------------------------------------------------------------
[PB Info 2022-Sep-21 16:18:02] || Program: Sorting Phase-II ||
[PB Info 2022-Sep-21 16:18:02] || Version: 4.0.0-1 ||
[PB Info 2022-Sep-21 16:18:02] || Start Time: Wed Sep 21 16:17:42 2022 ||
[PB Info 2022-Sep-21 16:18:02] || End Time: Wed Sep 21 16:18:02 2022 ||
[PB Info 2022-Sep-21 16:18:02] || Total Time: 20 seconds ||
[PB Info 2022-Sep-21 16:18:02] ------------------------------------------------------------------------------
[PB Info 2022-Sep-21 16:18:02] ------------------------------------------------------------------------------
[PB Info 2022-Sep-21 16:18:02] || Parabricks accelerated Genomics Pipeline ||
[PB Info 2022-Sep-21 16:18:02] || Version 4.0.0-1 ||
[PB Info 2022-Sep-21 16:18:02] || Marking Duplicates, BQSR ||
[PB Info 2022-Sep-21 16:18:02] ------------------------------------------------------------------------------
[PB Info 2022-Sep-21 16:18:02] progressMeter - Percentage
[PB Info 2022-Sep-21 16:18:12] 22.7 2.06 GB
[PB Info 2022-Sep-21 16:18:22] 44.1 2.16 GB
[PB Info 2022-Sep-21 16:18:32] 68.0 5.88 GB
[PB Info 2022-Sep-21 16:18:42] 96.0 0.20 GB
[PB Info 2022-Sep-21 16:18:52] 100.0 0.00 GB
[PB Info 2022-Sep-21 16:18:52] BQSR and writing final BAM: 50.370 seconds
[PB Info 2022-Sep-21 16:18:52] ------------------------------------------------------------------------------
[PB Info 2022-Sep-21 16:18:52] || Program: Marking Duplicates, BQSR ||
[PB Info 2022-Sep-21 16:18:52] || Version: 4.0.0-1 ||
[PB Info 2022-Sep-21 16:18:52] || Start Time: Wed Sep 21 16:18:02 2022 ||
[PB Info 2022-Sep-21 16:18:52] || End Time: Wed Sep 21 16:18:52 2022 ||
[PB Info 2022-Sep-21 16:18:52] || Total Time: 50 seconds ||
[PB Info 2022-Sep-21 16:18:52] ------------------------------------------------------------------------------

real 7m46.118s
user 107m4.313s
sys 2m16.551s

[user@cn4464 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226

```





