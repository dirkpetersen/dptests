

document.querySelector('title').textContent = 'fuseq-wes: discovering fusion genes from whole exome sequencing data in cancer patients';
**fuseq-wes: discovering fusion genes from whole exome sequencing data in cancer patients**


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



This tool is developed based on [FuSeq](https://hpc.nih.gov/apps/FuSeq.html), the method for detecting fusion genes from RNA-seq data. A subsampling study of the prostate data suggests that a coverage of at least 75x is necessary to achieve high accuracy.

### References:



Deng W, Murugan S, Lindberg J, Chellappa V, Shen X, Pawitan Y, Vu TN.   

*Fusion Gene Detection Using Whole-Exome Sequencing Data in Cancer Patients*    

[PubMed](https://pubmed.ncbi.nlm.nih.gov/35251131/)
[Front Genet](https://www.frontiersin.org/articles/10.3389/fgene.2022.820493/full), 2022, **13** 820493



Documentation
* [fuseq-wes GitHub Page](https://github.com/nghiavtr/FuSeq_WES)


Important Notes
* Module Name: fuseq-wes (see [the modules page](/apps/modules.html) for more information)
* Unusual environment variables set
	+ **FUSEQ\_WES** fuseq-wes installation directory
	+ **FUSEQ\_WES\_REF** fuseq-wes reference directory
	+ **FUSEQ\_WES\_TEST\_DATA** sample data for running fuseq-wes



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=4g --gres=lscratch:10** 
[user@cn3144 ~]$ **module load fuseq-wes** 
[+] Loading python 3.8  ...
[+] Loading gcc  9.2.0  ...
[-] Unloading gcc  9.2.0  ...
[+] Loading gcc  9.2.0  ...
[+] Loading openmpi 4.0.5  for GCC 9.2.0
[+] Loading ImageMagick  7.0.8  on cn4313
[+] Loading HDF5  1.10.4
[-] Unloading gcc  9.2.0  ...
[+] Loading gcc  9.2.0  ...
[+] Loading NetCDF 4.7.4_gcc9.2.0
[+] Loading pandoc  2.17.1.1  on cn4313
[+] Loading pcre2 10.21  ...
[+] Loading R 4.2.0
[+] Loading fuseq-wes  1.0.0

```

Create soft links to the sample read data:

```

[user@cn3144 ]$ **cp -r $FUSEQ\_WES\_TEST\_DATA/\* .** 
[user@cn3144 ]$ **bamfile="FuSeq\_WES\_testdata/test.bam"**
[user@cn3144 ]$ **ref\_json="$FUSEQ\_WES\_REF/UCSC\_hg19\_wes\_contigSize3000\_bigLen130000\_r100/\
UCSC\_hg19\_wes\_contigSize3000\_bigLen130000\_r100.json"**
[user@cn3144 ]$ **gtfSqlite="$FUSEQ\_WES\_REF/UCSC\_hg19\_wes\_contigSize3000\_bigLen130000\_r100/\
UCSC\_hg19\_wes\_contigSize3000\_bigLen130000\_r100.sqlite"**
[user@cn3144 ]$ **output\_dir="test\_out"**
[user@cn3144 ]$ **mkdir $output\_dir**

```

#extract mapped reads and split reads

```

[user@cn3144 ]$ **python3 $FUSEQ\_WES/FuSeq\_WES\_v1.0.0/fuseq\_wes.py \
 --bam $bamfile \
 --gtf $ref\_json \
 --mapq-filter \
 --outdir $output\_dir**

```

#process the reads

```


[user@cn3144 ]$ **fusiondbFn="$FUSEQ\_WES/FuSeq\_WES\_v1.0.0/Data/Mitelman\_fusiondb.RData"**
[user@cn3144 ]$ **paralogdb="$FUSEQ\_WES/FuSeq\_WES\_v1.0.0/Data/ensmbl\_paralogs\_grch37.RData"**
[user@cn3144 ]$ **Rscript $FUSEQ\_WES/FuSeq\_WES\_v1.0.0/process\_fuseq\_wes.R \
 in=$output\_dir \
 sqlite=$gtfSqlite \
 fusiondb=$fusiondbFn \
 paralogdb=$paralogdbFn \
 out=$output\_dir**

```

Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. fuseq-wes.sh). For example:



```

#! /bin/bash
module load fuseq-wes
set -e
cp -r $FUSEQ_WES_TEST_DATA/* .
bamfile="FuSeq_WES_testdata/test.bam"
ref_json="$FUSEQ_WES_REF/UCSC_hg19_wes_contigSize3000_bigLen130000_r100/\
UCSC_hg19_wes_contigSize3000_bigLen130000_r100.json"
gtfSqlite="$FUSEQ_WES_REF/UCSC_hg19_wes_contigSize3000_bigLen130000_r100/\
UCSC_hg19_wes_contigSize3000_bigLen130000_r100.sqlite"
output_dir="test_out"
mkdir -p $output_dir

python3 $FUSEQ_WES/FuSeq_WES_v1.0.0/fuseq_wes.py \
                   --bam $bamfile \
                   --gtf $ref_json \
                   --mapq-filter \
                   --outdir $output_dir


```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch -c 2 --mem=4g --time=8:00:00 fuseq-wes.sh
```

The master process submitting jobs should be run
either as a batch job or on an interactive node - not on the biowulf
login node.








