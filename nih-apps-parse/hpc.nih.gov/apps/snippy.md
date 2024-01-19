

document.querySelector('title').textContent = 'snippy: rapid haploid variant calling and core genome alignment ';
**snippy: rapid haploid variant calling and core genome alignment** 


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



snippy is a tool for rapid haploid variant calling and core genome alignment.




Documentation
* [snippy GitHub page](https://github.com/tseemann/snippy)
* [snippy tutorial](http://sepsis-omics.github.io/tutorials/modules/snippy/)


Important Notes
* Module Name: snippy (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set
	+ **SNIPPY\_HOME**  installation directory
	+ **SNIPPY\_BIN**       executable directory
	+ **SNIPPY\_SRC**       source code directory
	+ **SNIPPY\_DATA**  sample data directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=16g --gres=lscratch:10**
[user@cn3200 ~]$**module load snippy** 
[+] Loading snippy 4.4.1  ...

```

Here is how snippy can be run on a test example:

```

[user@cn3200 ~]$ **cp $SNIPPY\_DATA/\* .**
[user@cn3200 ~]$ **ln -s /fdb/igenomes/Homo\_sapiens/UCSC/hg19/hg19.fa**
[user@cn3200 ~]$ **snippy --cpus 16 --outdir mysnps --ref hg19.fa --R1 test\_R1.fastq.gz --R2 test\_R2.fastq.gz**
...
[11:14:32] Written by Torsten Seemann 
[11:14:32] Obtained from https://github.com/tseemann/snippy
[11:14:32] Detected operating system: linux
[11:14:32] Enabling bundled linux tools.
[11:14:32] Found bwa - /opt/conda/envs/snippy/bin/bwa
[11:14:32] Found samtools - /opt/conda/envs/snippy/bin/samtools
[11:14:32] Found tabix - /opt/conda/envs/snippy/bin/tabix
[11:14:32] Found bgzip - /opt/conda/envs/snippy/bin/bgzip
[11:14:32] Found snpEff - /opt/conda/envs/snippy/bin/snpEff
[11:14:32] Found java - /opt/conda/envs/snippy/bin/java
[11:14:32] Found gzip - /usr/bin/gzip
[11:14:32] Found parallel - /opt/conda/envs/snippy/bin/parallel
[11:14:32] Found freebayes - /opt/conda/envs/snippy/bin/freebayes
[11:14:32] Found freebayes-parallel - /opt/conda/envs/snippy/bin/freebayes-parallel
[11:14:32] Found fasta\_generate\_regions.py - /opt/conda/envs/snippy/bin/fasta\_generate\_regions.py
[11:14:32] Found vcfstreamsort - /opt/conda/envs/snippy/bin/vcfstreamsort
[11:14:32] Found vcfuniq - /opt/conda/envs/snippy/bin/vcfuniq
[11:14:32] Found vcffirstheader - /opt/conda/envs/snippy/bin/vcffirstheader
[11:14:32] Found vcf-consensus - /opt/conda/envs/snippy/bin/../binaries/noarch/vcf-consensus
[11:14:32] Found snippy-vcf\_to\_tab - /opt/conda/envs/snippy/bin/snippy-vcf\_to\_tab
[11:14:32] Found snippy-vcf\_report - /opt/conda/envs/snippy/bin/snippy-vcf\_report
[11:14:32] Found snippy-vcf\_filter - /opt/conda/envs/snippy/bin/snippy-vcf\_filter
[11:14:32] Checking version: samtools --version is >= 1.3 - ok, have 1.6
[11:14:32] Checking version: freebayes --version is >= 1.1 - ok, have 1.3
[11:14:32] Checking version: snpEff -version is >= 4.3 - ok, have 5.1
[11:14:32] Using reference: /gpfs/gsfs7/users/user/snippy/hg19.fa
[11:14:50] Treating reference as 'fasta' format.
[11:14:50] Will use 16 CPU cores.
[11:14:50] Using read file: /gpfs/gsfs7/users/user/snippy/test\_R1.fastq.gz
[11:14:50] Using read file: /gpfs/gsfs7/users/user/snippy/test\_R2.fastq.gz
[11:14:50] Output folder mysnps3 already exists.
[user@cn1113 snippy]$ snippy --cpus 16 --outdir mysnps4 --ref hg19.fa --R1 test\_R1.fastq.gz --R2 test\_R2.fastq.gz
[11:15:07] This is snippy 3.2-dev
[11:15:07] Written by Torsten Seemann 
[11:15:07] Obtained from https://github.com/tseemann/snippy
[11:15:07] Detected operating system: linux
[11:15:07] Enabling bundled linux tools.
[11:15:07] Found bwa - /opt/conda/envs/snippy/bin/bwa
[11:15:07] Found samtools - /opt/conda/envs/snippy/bin/samtools
[11:15:07] Found tabix - /opt/conda/envs/snippy/bin/tabix
[11:15:07] Found bgzip - /opt/conda/envs/snippy/bin/bgzip
[11:15:07] Found snpEff - /opt/conda/envs/snippy/bin/snpEff
[11:15:07] Found java - /opt/conda/envs/snippy/bin/java
[11:15:07] Found gzip - /usr/bin/gzip
[11:15:07] Found parallel - /opt/conda/envs/snippy/bin/parallel
[11:15:07] Found freebayes - /opt/conda/envs/snippy/bin/freebayes
[11:15:07] Found freebayes-parallel - /opt/conda/envs/snippy/bin/freebayes-parallel
[11:15:07] Found fasta\_generate\_regions.py - /opt/conda/envs/snippy/bin/fasta\_generate\_regions.py
[11:15:07] Found vcfstreamsort - /opt/conda/envs/snippy/bin/vcfstreamsort
[11:15:07] Found vcfuniq - /opt/conda/envs/snippy/bin/vcfuniq
[11:15:07] Found vcffirstheader - /opt/conda/envs/snippy/bin/vcffirstheader
[11:15:07] Found vcf-consensus - /opt/conda/envs/snippy/bin/../binaries/noarch/vcf-consensus
[11:15:07] Found snippy-vcf\_to\_tab - /opt/conda/envs/snippy/bin/snippy-vcf\_to\_tab
[11:15:07] Found snippy-vcf\_report - /opt/conda/envs/snippy/bin/snippy-vcf\_report
[11:15:07] Found snippy-vcf\_filter - /opt/conda/envs/snippy/bin/snippy-vcf\_filter
[11:15:07] Checking version: samtools --version is >= 1.3 - ok, have 1.6
[11:15:07] Checking version: freebayes --version is >= 1.1 - ok, have 1.3
[11:15:08] Checking version: snpEff -version is >= 4.3 - ok, have 5.1
[11:15:08] Using reference: /gpfs/gsfs7/users/user/snippy/hg19.fa
[11:15:19] Treating reference as 'fasta' format.
[11:15:19] Will use 16 CPU cores.
[11:15:19] Using read file: /gpfs/gsfs7/users/user/snippy/test\_R1.fastq.gz
[11:15:19] Using read file: /gpfs/gsfs7/users/user/snippy/test\_R2.fastq.gz
[11:15:19] Creating folder: mysnps4
[11:15:19] Changing working directory: mysnps4
[11:15:19] Creating reference folder: reference
[11:15:19] Extracting FASTA and GFF from reference.
[11:16:01] Wrote 25 sequences to ref.fa
[11:16:01] Wrote 0 features to ref.gff
[11:16:01] Freebayes will process 64 chunks of 49176391 bp, 16 chunks at a time.
[11:16:01] Using BAM RG (Read Group) ID: snps
[11:16:01] Running: samtools faidx reference/ref.fa 2>> snps.log
[11:16:17] Running: bwa index reference/ref.fa 2>> snps.log

[12:01:52] Running: mkdir reference/genomes && cp -f reference/ref.fa reference/genomes/ref.fa 2>> snps.log
[12:01:54] Running: mkdir reference/ref && bgzip -c reference/ref.gff > reference/ref/genes.gff.gz 2>> snps.log
[12:01:54] Running: (bwa mem -v 2 -M -R '@RG\tID:snps\tSM:snps' -t 16 reference/ref.fa /gpfs/gsfs7/users/user/snippy/test\_R1.fastq.gz /gpfs/gsfs7/users/user/snippy/test\_R2.fastq.gz | samtools view -u -T reference/ref.fa -F 3844 -q 60 | samtools sort --reference reference/ref.fa > snps.bam) 2>> snps.log
[12:01:57] Running: samtools index snps.bam 2>> snps.log
[12:01:58] Running: samtools depth -aa -q 20 snps.bam | bgzip > snps.depth.gz 2>> snps.log
[12:07:08] Running: tabix -s 1 -b 2 -e 2 snps.depth.gz 2>> snps.log
[12:11:06] Running: fasta\_generate\_regions.py reference/ref.fa.fai 49176391 > reference/ref.txt 2>> snps.log
[12:11:07] Running: freebayes-parallel reference/ref.txt 16 -p 1 -q 20 -m 60 --min-coverage 10 -V -f reference/ref.fa snps.bam > snps.raw.vcf 2>> snps.log
[12:11:25] Running: /opt/conda/envs/snippy/bin/snippy-vcf\_filter --minqual 10 --mincov 10 --minfrac 0.9 snps.raw.vcf > snps.filt.vcf 2>> snps.log
[12:11:25] Running: cp snps.filt.vcf snps.vcf 2>> snps.log
[12:11:25] Running: bgzip -c snps.vcf > snps.vcf.gz 2>> snps.log
[12:11:25] Running: tabix -p vcf snps.vcf.gz 2>> snps.log
[12:11:25] Running: /opt/conda/envs/snippy/bin/snippy-vcf\_to\_tab --gff reference/ref.gff --ref reference/ref.fa --vcf snps.vcf > snps.tab 2>> snps.log
[12:11:49] Running: vcf-consensus snps.vcf.gz < reference/ref.fa > snps.consensus.fa 2>> snps.log
[12:13:43] Running: /opt/conda/envs/snippy/bin/snippy-vcf\_filter --subs snps.filt.vcf > snps.filt.subs.vcf 2>> snps.log
[12:13:43] Running: bgzip -c snps.filt.subs.vcf > snps.filt.subs.vcf.gz 2>> snps.log
[12:13:43] Running: tabix -p vcf snps.filt.subs.vcf.gz 2>> snps.log
[12:13:43] Running: vcf-consensus snps.filt.subs.vcf.gz < reference/ref.fa > snps.consensus.subs.fa 2>> snps.log
[12:15:40] Generating aligned/masked FASTA relative to reference: snps.aligned.fa

...

[12:46:32] Creating extra output files: BED GFF CSV TXT HTML
[12:46:32] Identified 3 variants.
[12:46:32] Result folder: mysnps4
[12:46:32] Result files:
[12:46:32] \* mysnps4/snps.aligned.fa
[12:46:32] \* mysnps4/snps.bam
[12:46:32] \* mysnps4/snps.bam.bai
[12:46:32] \* mysnps4/snps.bed
[12:46:32] \* mysnps4/snps.consensus.fa
[12:46:32] \* mysnps4/snps.consensus.subs.fa
[12:46:32] \* mysnps4/snps.csv
[12:46:32] \* mysnps4/snps.depth.gz
[12:46:32] \* mysnps4/snps.depth.gz.tbi
[12:46:32] \* mysnps4/snps.filt.subs.vcf
[12:46:32] \* mysnps4/snps.filt.subs.vcf.gz
[12:46:32] \* mysnps4/snps.filt.subs.vcf.gz.tbi
[12:46:32] \* mysnps4/snps.filt.vcf
[12:46:32] \* mysnps4/snps.gff
[12:46:32] \* mysnps4/snps.html
[12:46:32] \* mysnps4/snps.log
[12:46:32] \* mysnps4/snps.raw.vcf
[12:46:32] \* mysnps4/snps.tab
[12:46:32] \* mysnps4/snps.txt
[12:46:32] \* mysnps4/snps.vcf
[12:46:32] \* mysnps4/snps.vcf.gz
[12:46:32] \* mysnps4/snps.vcf.gz.tbi
[12:46:32] Walltime used: 1 hours, 31 minutes, 25 seconds
[12:46:32] Found a bug? Post it at https://github.com/tseemann/snippy/issues
[12:46:32] Done.

```

An output folder **mysnps**  is created:

```

[user@cn3200 ~]$ **tree mysnps**
mysnps
├ reference
│   ├ genomes
│   │   └ ref.fa
│   ├ ref
│   │   ├ genes.gff.gz
│   │   └ snpEffectPredictor.bin
│   ├ ref.fa
│   ├ ref.fa.amb
│   ├ ref.fa.ann
│   ├ ref.fa.bwt
│   ├ ref.fa.fai
│   ├ ref.fa.pac
│   ├ ref.fa.sa
│   ├ ref.gff
│   ├ ref.txt
│   └ snpeff.config
├ ref.fa -> reference/ref.fa
├ ref.fa.fai -> reference/ref.fa.fai
├ snps.aligned.fa
├ snps.bam
├ snps.bam.bai
├ snps.bed
├ snps.consensus.fa
├ snps.consensus.subs.fa
├ snps.csv
├ snps.filt.vcf
├ snps.gff
├ snps.html
├ snps.log
├ snps.raw.vcf
├ snps.subs.vcf
├ snps.tab
├ snps.txt
├ snps.vcf
├ snps.vcf.gz
└ snps.vcf.gz.csi

3 directories, 33 filesa

```

End the interactive session:

```

[user@cn3200 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```





