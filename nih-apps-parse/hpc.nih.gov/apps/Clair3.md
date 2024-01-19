

document.querySelector('title').textContent = 'Clair3: Integrating pileup and full-alignment for high-performance long-read variant calling';
**Clair3: Integrating pileup and full-alignment for high-performance long-read variant calling**


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



Clair3 is a small variant caller for Illumina, PacBio and ONT long reads. Compare to PEPPER (r0.4), Clair3 (v0.1) shows a better SNP F1-score with ≤30-fold of ONT data (precisionFDA Truth Challenge V2), and a better Indel F1-score, while runs generally four times faster.



Documentation
* [Clair3 Github / tutorial page](https://github.com/HKU-BAL/Clair3)
* [Training data and trained models](https://github.com/HKU-BAL/Clair3/blob/main/docs/training_data.md)
* [Demos for processing Illumina, PacBio and ONT data](https://github.com/HKU-BAL/Clair3/tree/main/docs/quick_demo)


Important Notes
* Module Name: Clair3 (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Multithreaded
* Can be run on both CPU and GPU nodes
* Unusual environment variables set
	+ **CLAIR3\_HOME**  installation directory
	+ **CLAIR3\_BIN**       executable directory
	+ **CLAIR3\_MODELS**    trained models folder
	+ **CLAIR3\_DATA**    sample data folder



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. 
  
Sample session on a GPU node:



```

[user@biowulf ~]$ **sinteractive --gres=gpu:p100:1,lscratch:10 --mem=16g -c4**
[user@cn2379 ~]$ **module load clair3**
[user@cn2379 ~]$ **cp -r $CLAIR3\_DATA/\* .**

[user@cn2379 ~]$ **THREADS=4**
[user@cn2379 ~]$ **OUTPUT\_VCF\_FILE\_PATH=merge\_output.vcf.gz**

```

 **Processing sample Illumina data**  


```
 
[user@cn2379 ~]$ **PLATFORM='ilmn'** 
[user@cn2379 ~]$ **INPUT\_DIR="Illumina"
[user@cn2379 ~]$ **cp -r $CLAIR3\_DATA/${INPUT\_DIR} .**
[user@cn2379 ~]$ **REF="GRCh38\_chr20.fa"**
[user@cn2379 ~]$ **BAM="HG003\_chr20\_demo.bam"**
[user@cn2379 ~]$ **BASELINE\_VCF\_FILE\_PATH="HG003\_GRCh38\_chr20\_v4.2.1\_benchmark.vcf.gz"**
[user@cn2379 ~]$ **BASELINE\_BED\_FILE\_PATH="HG003\_GRCh38\_chr20\_v4.2.1\_benchmark\_noinconsistent.bed"**
[user@cn2379 ~]$ **clair3 \
 --bam\_fn=${INPUT\_DIR}/${BAM} \
 --ref\_fn=${INPUT\_DIR}/${REF} \
 --model\_path=${CLAIR3\_MODELS}/${PLATFORM} \
 --threads=${THREADS} \
 --platform=${PLATFORM} \
 --output=./ \
 --bed\_fn=${INPUT\_DIR}/${BASELINE\_BED\_FILE\_PATH}** 
...**
```

 **Processing sample PacBio Hifi data**  


```

[user@cn2379 ~]$ **PLATFORM='hifi'** 
[user@cn2379 ~]$ **INPUT\_DIR="PacBio"
[user@cn2379 ~]$ **cp -r $CLAIR3\_DATA/${INPUT\_DIR} .**
[user@cn2379 ~]$ **REF="GRCh38\_no\_alt\_chr20.fa"**
[user@cn2379 ~]$ **BAM="HG003\_chr20\_demo.bam"**
[user@cn2379 ~]$ **BASELINE\_VCF\_FILE\_PATH="HG003\_GRCh38\_chr20\_v4.2.1\_benchmark.vcf.gz"**
[user@cn2379 ~]$ **BASELINE\_BED\_FILE\_PATH="HG003\_GRCh38\_chr20\_v4.2.1\_benchmark\_noinconsistent.bed"**
[user@cn2379 ~]$ **clair3 \
 --bam\_fn=${INPUT\_DIR}/${BAM} \
 --ref\_fn=${INPUT\_DIR}/${REF} \
 --threads=${THREADS} \
 --platform=${PLATFORM} \
 --model\_path=${CLAIR3\_MODELS}/${PLATFORM} \
 --output=./ \
 --bed\_fn=${INPUT\_DIR}/${BASELINE\_BED\_FILE\_PATH}** 
...**
```

 **Processing sample ONT data**  


```

[user@cn2379 ~]$ **PLATFORM='ont'** 
[user@cn2379 ~]$ **INPUT\_DIR="ONT"
[user@cn2379 ~]$ **cp -r $CLAIR3\_DATA/${INPUT\_DIR} .**
[user@cn2379 ~]$ **REF="GRCh38\_no\_alt\_chr20.fa"**
[user@cn2379 ~]$ **BAM="HG003\_chr20\_demo.bam"**
[user@cn2379 ~]$ **BASELINE\_VCF\_FILE\_PATH="HG003\_GRCh38\_chr20\_v4.2.1\_benchmark.vcf.gz"**
[user@cn2379 ~]$ **BASELINE\_BED\_FILE\_PATH="HG003\_GRCh38\_chr20\_v4.2.1\_benchmark\_noinconsistent.bed"**
[user@cn2379 ~]$ **clair3 \
 --bam\_fn=${INPUT\_DIR}/${BAM} \
 --ref\_fn=${INPUT\_DIR}/${REF} \
 --threads=${THREADS} \
 --platform=${PLATFORM} \ --model\_path=${CLAIR3\_MODELS}/${PLATFORM} \
 --output=./ \
 --vcf\_fn=${INPUT\_DIR}/${BASELINE\_VCF\_FILE\_PATH}** 
...**
```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. clair3.sh). For example:



```

#!/bin/bash
set -e
module load Clair3 
...               
clair3 \ 
  --bam_fn=${INPUT_DIR}/${BAM} \
  --ref_fn=${INPUT_DIR}/${REF} \                                                                      --threads=${THREADS} \
  --platform=${PLATFORM} \                                                                            --model_path=${CLAIR3_MODELS}/${PLATFORM} \
  --output=./ \
  --vcf_fn=${INPUT_DIR}/${BASELINE_VCF_FILE_PATH}

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
**sbatch clair3.sh**
```







