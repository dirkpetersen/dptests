

document.querySelector('title').textContent = 'crossmap on Biowulf';
crossmap on Biowulf


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



CrossMap is a program for convenient conversion of genome coordinates between different assemblies (e.g. mm9->mm10). It can convert SAM, BAM, bed, GTF, GFF, wig/bigWig, and VCF files.



### References:


* Hao Zhao *et al.* *CrossMap: a versatile tool for coordinate 
 conversion between genome assemblies*. Bioinformatics 2014(30): 1006-1007.
 [PubMed](http://www.ncbi.nlm.nih.gov/pubmed/24351709)  | 
 [PMC](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3967108/) | 
 [Journal](http://bioinformatics.oxfordjournals.org/content/30/7/1006)


Documentation
* [crossmap Main Site](http://crossmap.sourceforge.net/)


Important Notes
* Module Name: crossmap (see [the modules page](/apps/modules.html) for more information)
* Singlethreaded app
* Example files in /usr/local/apps/crossmap/TEST\_DATA



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session (user input in **bold**):



```

[user@biowulf ~]$ **sinteractive -c2 --mem=4g --gres=lscratch:10**
salloc.exe: Pending job allocation 11342506
salloc.exe: job 11342506 queued and waiting for resources
salloc.exe: job 11342506 has been allocated resources
salloc.exe: Granted job allocation 11342506
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn0865 are ready for job
srun: error: x11: no local DISPLAY defined, skipping
error: unable to open file /tmp/slurm-spank-x11.11342506.0
slurmstepd: error: x11: unable to read DISPLAY value

[user@cn0865 ~]$ **cd /lscratch/$SLURM\_JOB\_ID**

[user@cn0865 11342506]$ **module load crossmap**
[+] Loading crossmap  0.5.2  on cn0865
[+] Loading singularity  3.7.2  on cn0865

[user@cn0865 11342506]$ **cp $CMAP\_DATA/\* .**

[user@cn0865 11342506]$ **crossmap -h**
Program: CrossMap (v0.2.8)

Description:
  CrossMap is a program for convenient conversion of genome coordinates and genome
  annotation files between assemblies (eg. lift from human hg18 to hg19 or vice
  versa). It supports file in BAM, SAM, BED, Wiggle, BigWig, GFF, GTF and VCF
  format.

Usage: CrossMap.py  [options]

  bam   convert alignment file in BAM or SAM format.
  bed   convert genome cooridnate or annotation file in BED or BED-like format.
  bigwig        convert genome coordinate file in BigWig format.
  gff   convert genome cooridnate or annotation file in GFF or GTF format.
  vcf   convert genome coordinate file in VCF format.
  wig   convert genome coordinate file in Wiggle, or bedGraph format.


[user@cn0865 11342506]$ **crossmap bed hg18ToHg19.over.chain test\_input > test\_output**
@ 2021-03-25 12:17:39: Read chain_file:  hg18ToHg19.over.chain

[user@cn0865 11342506]$ **head -n2 test\_output**
chr1    142614848       142617697       ->      chr1    143903503       143906352
chr1    142617697       142623312       ->      chr1    143906355       143911970

[user@cn0865 11342506]$ **diff --ignore-all-space expected\_output test\_output**

[user@cn0865 11342506]$ **exit**
exit
salloc.exe: Relinquishing job allocation 11342506

[user@biowulf ~]$ 

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. crossmap.sh). For example:



```

#!/bin/bash
function fail() {
    echo "$@" >&2
    exit 1
}

module load crossmap || fail "could not load crossmap module"
if [[ ! -f hg19ToHg38.over.chain.gz ]]; then
    wget http://hgdownload.soe.ucsc.edu/goldenPath/mm9/liftOver/hg19ToHg38.over.chain.gz
fi
crossmap bam hg19ToHg38.over.chain.gz hg19_example.bam out

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch crossmap.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. crossmap.swarm). For example:



```

crossmap bam hg19ToHg38.over.chain.gz sample1.bam sample1_hg38.bam
crossmap bam hg19ToHg38.over.chain.gz sample2.bam sample2_hg38.bam
crossmap bam hg19ToHg38.over.chain.gz sample3.bam sample3_hg38.bam

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f crossmap.swarm --module crossmap
```

where


|  |  |
| --- | --- |
| --module crossmap Loads the crossmap module for each subjob in the swarm 
 | |








