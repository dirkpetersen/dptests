

document.querySelector('title').textContent = 'smoove: structural variant calling and genotyping with existing tools, but, smoothly ';
**smoove: structural variant calling and genotyping with existing tools, but, smoothly** 


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



smoove simplifies and speeds calling and genotyping SVs for short reads. 
It also improves specificity by removing many spurious alignment signals 
that are indicative of low-level noise and often contribute to spurious calls.



Documentation
* [smoove GitHub page](https://github.com/brentp/smoove)


Important Notes
* Module Name: smoove (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set
	+ **SMOOVE\_HOME**  installation directory
	+ **SMOOVE\_BIN**       executable directory
	+ **SMOOVE\_DATA**  test data directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=24g -c12 --gres=lscratch:50**
[user@cn3200 ~]$**module load smoove** 
[+] Loading singularity  3.5.3  on cn3112
[+] Loading smoove 0.2.7  ...

[user@cn3200 ~]$ **mkdir out\_dir**
[user@cn3200 ~]$ **export TMPDIR=`pwd`**

```

Copy sample data to the current folder:

```

[user@cn3200 ~]$ **cp -P $SMOOVE\_DATA/\* .**

```

Trim reads:

```

[user@cn3200 ~]$ **trimmomatic PE -threads 12 ./sample.R1.fastq ./sample.R2.fastq ./sample.R1.paired.fastq ./sample.R1.unpaired.fastq ./sample.R2.paired.fastq ./sample.R2.unpaired.fastq ILLUMINACLIP:/usr/local/apps/trimmomatic/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36**
Using PrefixPair: 'TACACTCTTTCCCTACACGACGCTCTTCCGATCT' and 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT'
ILLUMINACLIP: Using 1 prefix pairs, 0 forward/reverse sequences, 0 forward only sequences, 0 reverse only sequences
Quality encoding detected as phred33
Input Read Pairs: 500000 Both Surviving: 488288 (97.66%) Forward Only Surviving: 8039 (1.61%) Reverse Only Surviving: 2828 (0.57%) Dropped: 845 (0.17%)
TrimmomaticPE: Completed successfully

```

From the resulting paired reads, produce BAM files sorted by position:

```

[user@cn3200 ~]$ **ln -s /fdb/igenomes/Homo\_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa hg38.fa**
[user@cn3200 ~]$ **ln -s /fdb/igenomes/Homo\_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa.fai hg38.fa.fai**
[user@cn3200 ~]$ **bwa index hg38.fa**
[bwa_index] Pack FASTA... 30.16 sec
[bwa_index] Construct BWT for the packed sequence...
[BWTIncCreate] textLength=6199845082, availableWord=448243540
[BWTIncConstructFromPacked] 10 iterations done. 99999994 characters processed.
[BWTIncConstructFromPacked] 20 iterations done. 199999994 characters processed.
...
[BWTIncConstructFromPacked] 670 iterations done. 6160559482 characters processed.
[BWTIncConstructFromPacked] 680 iterations done. 6184133946 characters processed.
[bwt_gen] Finished constructing BWT in 688 iterations.
[bwa_index] 3436.59 seconds elapse.
[bwa_index] Update BWT... 20.91 sec
[bwa_index] Pack forward-only FASTA... 19.82 sec
[bwa_index] Construct SA from BWT and Occ...
1085.34 sec
[main] Version: 0.7.17-r1188
[main] CMD: /opt/conda/envs/smoove-env/bin/bwa index hg38.fa
[main] Real time: 4608.042 sec; CPU: 4592.900 sec
[user@cn3200 ~]$ **bwa mem -t 12 hg38.fa sample.R1.paired.fastq sample.R2.paired.fastq | samtools sort -o out\_dir/sample\_sorted.bam**
[M::bwa_idx_load_from_disk] read 0 ALT contigs
[M::process] read 816658 sequences (120000292 bp)...
[M::process] read 159918 sequences (23530318 bp)...
 [M::mem_pestat] # candidate unique pairs for (FF, FR, RF, RR): (2852, 227554, 372, 3105)
[M::mem_pestat] analyzing insert size distribution for orientation FF...
[M::mem_pestat] (25, 50, 75) percentile: (632, 1245, 2642)
[M::mem_pestat] low and high boundaries for computing mean and std.dev: (1, 6662)
[M::mem_pestat] mean and std.dev: (1672.47, 1532.16)
[M::mem_pestat] low and high boundaries for proper pairs: (1, 8672)
...
[M::mem_pestat] skip orientation FF
[M::mem_pestat] skip orientation RF
[M::mem_pestat] skip orientation RR
[M::mem_process_seqs] Processed 159918 reads in 337.782 CPU sec, 85.045 real sec
[main] Version: 0.7.17-r1188
[main] CMD: bwa mem -t 12 hg38.fa sample.R1.paired.fastq sample.R2.paired.fastq
[main] Real time: 547.623 sec; CPU: 2046.621 sec

```

As needed, add a missing read group (RG) tag to the BAM file:

```

[user@cn3200 ~]$ **picard AddOrReplaceReadGroups I=out\_dir/sample\_sorted.bam O=out\_dir/sample\_sorted\_rg.bam RGID=4 \
RGLB=lib1 \
RGPL=ILLUMINA \
RGPU=unit1 \
RGSM=20**
INFO    2020-06-01 12:51:12     AddOrReplaceReadGroups

********** NOTE: Picard's command line syntax is changing.
**********
********** For more information, please see:
********** https://github.com/broadinstitute/picard/wiki/Command-Line-Syntax-Transition-For-Users-(Pre-Transition)
**********
********** The command line looks like this in the new syntax:
**********
**********    AddOrReplaceReadGroups -I sample_sorted.bam -O sample_sorted_p_rg.bam -RGID 4 -RGLB lib1 -RGPL ILLUMINA -RGPU unit1 -RGSM 20
**********


12:51:12.787 INFO  NativeLibraryLoader - Loading libgkl_compression.so from jar:file:/usr/local/apps/picard/2.22.2/picard.jar!/com/intel/gkl/native/libgkl_compression.so
[Mon Jun 01 12:51:12 EDT 2020] AddOrReplaceReadGroups INPUT=sample_sorted.bam OUTPUT=sample_sorted_p_rg.bam RGID=4 RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM=20    VERBOSITY=INFO QUIET=false VALIDATION_STRINGENCY=STRICT COMPRESSION_LEVEL=5 MAX_RECORDS_IN_RAM=500000 CREATE_INDEX=false CREATE_MD5_FILE=false GA4GH_CLIENT_SECRETS=client_secrets.json USE_JDK_DEFLATER=false USE_JDK_INFLATER=false
[Mon Jun 01 12:51:12 EDT 2020] Executing as user@cn3133 on Linux 3.10.0-862.14.4.el7.x86_64 amd64; OpenJDK 64-Bit Server VM 1.8.0_181-b13; Deflater: Intel; Inflater: Intel; Provider GCS is not available; Picard version: 2.22.2
INFO    2020-06-01 12:51:12     AddOrReplaceReadGroups  Created read-group ID=4 PL=ILLUMINA LB=lib1 SM=20

INFO    2020-06-01 12:51:31     AddOrReplaceReadGroups  Processed     1,000,000 records.  Elapsed time: 00:00:18s.  Time for last 1,000,000:   18s.  Last read position: chrX:37,671,305
[Mon Jun 01 12:51:32 EDT 2020] picard.sam.AddOrReplaceReadGroups done. Elapsed time: 0.33 minutes.
Runtime.totalMemory()=2058354688

```

Index the resulting BAM file:

```

[user@cn3200 ~]$ **samtools index out\_dir/sample\_sorted\_rg.bam**

```

Run smoove on the BAM file:

```

[user@cn3200 ~]$ **smoove -h**
smoove version: 0.2.7

smoove calls several programs. Those with 'Y' are found on your $PATH. Only those with '*' are required.

 *[Y] bgzip [ sort   -> (compress) ->   index ]
 *[Y] gsort [(sort)  ->  compress   ->  index ]
 *[Y] tabix [ sort   ->  compress   -> (index)]
 *[Y] lumpy
 *[Y] lumpy_filter
 *[Y] samtools
 *[Y] svtyper
 *[Y] mosdepth [extra filtering of split and discordant files for better scaling]

  [Y] duphold [(optional) annotate calls with depth changes]
  [Y] svtools [only needed for large cohorts].

Available sub-commands are below. Each can be run with -h for additional help.

call     : call lumpy (and optionally svtyper)
merge    : merge and sort (using svtools) calls from multiple samples
genotype : parallelize svtyper on an input VCF
paste    : square final calls from multiple samples (each with same number of variants)
annotate : annotate a VCF with gene and quality of SV call
hipstr   : run hipSTR in parallel
cnvnator : run cnvnator in parallel
duphold  : run duphold in parallel (this can be done by adding a flag to call or genotype)

[user@cn3200 ~]$ **smoove call -o out\_dir --processes 2 --fasta hg38.fa --name sample\_sorted\_bam --excludechroms '~^GL,~^HLA,~\_random,~^chrUn,~alt,~decoy' out\_dir/sample\_sorted\_rg.bam**
[smoove] 2022/08/18 12:20:45 starting with version 0.2.7
[smoove] 2022/08/18 12:20:45 calculating bam stats for 1 bams
...
[smoove]:2022/08/18 12:20:48 finished process: lumpy-filter (set -eu; lumpy_filter -f hg38.fa sample_sorted.bam out_dir/sample_sorted.split.bam.tmp.bam out_dir/s) in user-time:3.810234s system-time:662.965ms
[smoove] 2022/08/18 12:20:48 done calculating bam stats
bash: warning: setlocale: LC_ALL: cannot change locale (en_US.UTF-8)
bash: warning: setlocale: LC_ALL: cannot change locale (en_US.UTF-8)
[smoove] 2022/08/18 12:21:24 removed 48530 alignments out of 128238 (37.84%) with low mapq, depth > 1000, or from excluded chroms from sample_sorted.split.bam in 36 seconds
[smoove] 2022/08/18 12:21:24 removed 5716 alignments out of 128238 (4.46%) that were bad interchromosomals or flanked-splitters from sample_sorted.split.bam
[smoove] 2022/08/18 12:21:25 removed 63824 alignments out of 207734 (30.72%) with low mapq, depth > 1000, or from excluded chroms from sample_sorted.disc.bam in 37 seconds
[smoove] 2022/08/18 12:21:25 removed 8109 alignments out of 207734 (3.90%) that were bad interchromosomals or flanked-splitters from sample_sorted.disc.bam
[smoove] 2022/08/18 12:21:27 kept 10426 putative orphans
[smoove] 2022/08/18 12:21:27 removed 37923 discordant orphans in 1 seconds
[smoove] 2022/08/18 12:21:28 removed 96013 singletons and isolated interchromosomals of 135801 reads (70.70%) from sample_sorted.disc.bam in 3 seconds
[smoove] 2022/08/18 12:21:28 39788 reads (19.15%) of the original 207734 remain from sample_sorted.disc.bam
[smoove] 2022/08/18 12:21:29 kept 6816 putative orphans
[smoove] 2022/08/18 12:21:29 removed 18823 split orphans in 4 seconds
[smoove] 2022/08/18 12:21:29 removed 46530 singletons of 73992 reads (62.89%) from sample_sorted.split.bam in 5 seconds
[smoove] 2022/08/18 12:21:29 27462 reads (21.41%) of the original 128238 remain from sample_sorted.split.bam
[smoove] 2022/08/18 12:21:29 starting lumpy
[smoove] 2022/08/18 12:21:29 wrote lumpy command to out_dir/sample_sorted.bam-lumpy-cmd.sh
[smoove] 2022/08/18 12:21:29 bash: warning: setlocale: LC_ALL: cannot change locale (en_US.UTF-8)
[smoove] 2022/08/18 12:21:29 chr1       1000000
chr1    2000000
[smoove] 2022/08/18 12:21:30 chr10      1000000
[smoove] 2022/08/18 12:21:30 chr11      1000000
[smoove] 2022/08/18 12:21:30 chr12      1000000
[smoove] 2022/08/18 12:21:30 chr13      1000000
chr13   2000000
[smoove] 2022/08/18 12:21:30
chr13   4000000
chr13   8000000
chr13   16000000
chr13   32000000
[smoove] 2022/08/18 12:21:30 chr14      1000000
[smoove] 2022/08/18 12:21:30 chr14      2000000
chr14   4000000
chr14   8000000
chr14   16000000
chr14   32000000
[smoove] 2022/08/18 12:21:30 chr15      1000000
chr15   2000000
[smoove] 2022/08/18 12:21:30 chr15      4000000
chr15   8000000
chr15   16000000
chr15
[smoove] 2022/08/18 12:21:30    32000000
[smoove] 2022/08/18 12:21:30 chr16      1000000
[smoove] 2022/08/18 12:21:30 chr17      1000000
[smoove] 2022/08/18 12:21:31 chr18      1000000
[smoove] 2022/08/18 12:21:31
[smoove] 2022/08/18 12:21:31 chr19      1000000
[smoove] 2022/08/18 12:21:31 chr2       1000000
[smoove] 2022/08/18 12:21:31 chr20      1000000
[smoove] 2022/08/18 12:21:31 chr21      1000000
chr21
[smoove] 2022/08/18 12:21:31    2000000
chr21   4000000
chr21   8000000
[smoove] 2022/08/18 12:21:31 chr22      1000000
chr22
[smoove] 2022/08/18 12:21:31    2000000
chr22   4000000
chr22   8000000
chr22   16000000
[smoove] 2022/08/18 12:21:31 chr3       1000000
[smoove] 2022/08/18 12:21:31
[smoove] 2022/08/18 12:21:31 chr4
[smoove] 2022/08/18 12:21:31    1000000
[smoove] 2022/08/18 12:21:32 chr5       1000000
[smoove] 2022/08/18 12:21:32 chr6       1000000
[smoove] 2022/08/18 12:21:32 chr7       1000000
[smoove] 2022/08/18 12:21:32 chr8       1000000
[smoove] 2022/08/18 12:21:32 chr8       4000000
[smoove] 2022/08/18 12:21:32 chr9       1000000
[smoove] 2022/08/18 12:21:33 chrM       1000000
[smoove] 2022/08/18 12:21:33 chrX       1000000
chrX
[smoove] 2022/08/18 12:21:33    2000000
chrX    4000000
[smoove] 2022/08/18 12:21:33 chrY       1000000
[smoove] 2022/08/18 12:21:33 chrY       2000000
chrY    4000000
chrY    8000000
chrY    16000000
[smoove] 2022/08/18 12:21:33 chrY       64000000
[smoove] 2022/08/18 12:21:36 wrote to out_dir/sample_sorted_bam-smoove.vcf.gz
[user@cn3200 ~]$ **smoove merge --fasta hg38.fa -o out\_dir -name sample\_sorted\_svs\_bam out\_dir/sample\_sorted\_bam-smoove.vcf.gz**
[smoove] 2022/08/18 13:16:27 starting with version 0.2.7
[smoove] 2022/08/18 13:16:27 merging 1 files
[smoove] 2022/08/18 13:16:27 finished sorting 1 files; merge starting.
[smoove] 2022/08/18 13:16:29 bash: warning: setlocale: LC_ALL: cannot change locale (en_US.UTF-8)
[smoove] 2022/08/18 13:16:30 wrote sites file to out_dir/sample_sorted_svs_bam.sites.vcf.gz
[smoove] 2022/08/18 13:16:30 wrote html file of disc, split counts to out_dir/sample_sorted_svs_bam.smoove-counts.html
[user@cn3200 ~]$ **smoove genotype -x --fasta hg38.fa -o out\_dir --name sample\_sorted\_svs\_bam-g --vcf out\_dir/sample\_sorted\_svs\_bam.sites.vcf.gz out\_dir/sample\_sorted\_rg.bam**
[smoove] 2022/08/18 13:28:48 starting with version 0.2.7
[smoove] 2022/08/18 13:28:48 writing sorted, indexed file to out_dir/sample_sorted_svs_bam-g-smoove.genotyped.vcf.gz
[smoove] 2022/08/18 13:28:48 bash: warning: setlocale: LC_ALL: cannot change locale (en_US.UTF-8)
[smoove] 2022/08/18 13:28:48 > gsort version 0.0.6
[smoove] 2022/08/18 13:28:52 wrote sorted, indexed file to out_dir/sample_sorted_svs_bam-g-smoove.genotyped.vcf.gz

```

End the interactive session:

```

[user@cn3200 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```

0.2.7




