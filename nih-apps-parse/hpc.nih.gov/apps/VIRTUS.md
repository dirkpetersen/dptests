

document.querySelector('title').textContent = 'VIRTUS on Biowulf';
VIRTUS on Biowulf


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



VIRTUS is a bioinformatics pipeline for viral transcriptome detection and quantification considering splicing.


The HPC version do not support docker, but use singularity instead.
Please do not install cwltool(3.1.20211014180718), since it can't work well with singularity CE.







Features
* VIRTUS is the first tool to detect viral transcripts considering their splicing event rather than the viral genome copy number. 
 * VIRTUS can be applied to both bulk RNAseq and single-cell RNAseq.
 * The virus reference covers 762 viruses including SARS-CoV-2 (cause of COVID-19).
 * The workflow is implemented by Common Workflow Language and Rabix. 
 * You can specify each parameter individually or give yaml or json file which describes all the parameter information.




### References:


* Yasumizu, Yoshiaki, Atsushi Hara, Shimon Sakaguchi, and Naganari Ohkura. */VIRTUS: a pipeline for comprehensive virus analysis from conventional RNA-seq data*Bioinformatics. 2020;btaa859. doi:10.1093/bioinformatics/btaa859
 [PubMed](https://www.ncbi.nlm.nih.gov/pubmed/33017003) | 
 [Journal](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btaa859/5918022)


Documentation
* VIRTUS Main Site:[Main Site](https://github.com/yyoshiaki/VIRTUS)
* VIRTUS2 Main Site:[Main Site](https://github.com/yyoshiaki/VIRTUS2)


Important Notes
* Module Name: VIRTUS (see [the modules page](/apps/modules.html) for more information)
 * VIRTUS is easy to run with wrapper:
 
```

	VIRTUS_wrapper.py -h
```
* VIRTUS workflow can also run with cwltool:
 
```

        cwltool --singularity --tmp-outdir-prefix=/lscratch/$SLURM_JOB_ID/ \
	--tmpdir-prefix=/lscratch/$SLURM_JOB_ID/ \
	$VIRTUS_WORKFLOW/VIRTUS.PE.cwl
```
* VIRTUS workflow has to be run with "--singularity" option. Using "--singularity" will download several containers to your local directory. 
 * You need to locate lscratch, because by default cwltools can easily fill up /tmp.
* Environment variables set 
	+ $VIRTUS\_WORKFLOW
	+ $VIRTUS\_INDEX* Example files in $VIRTUS\_TEST\_DATA
* VIRTUS/1 and VIRTUS/2 are running differently: see examples bellow.



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive --gres=lscratch:500 --cpus-per-task=40 --mem=40G**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load VIRTUS/1.2.1**
[user@cn3144 ~]$ **mkdir -p /data/$USER/VIRTUS; cd /data/$USER/VIRTUS**
[user@cn3144 VIRTUS]$ **cp "$VIRTUS\_TEST\_DATA"/\* .**
[user@cn3144 VIRTUS]$ **VIRTUS\_wrapper.py -h**
usage: VIRTUS_wrapper.py [-h] [--VIRTUSDir VIRTUSDIR] --genomeDir_human
                         GENOMEDIR_HUMAN --genomeDir_virus GENOMEDIR_VIRUS
                         --salmon_index_human SALMON_INDEX_HUMAN
                         [--salmon_quantdir_human SALMON_QUANTDIR_HUMAN]
                         [--outFileNamePrefix_human OUTFILENAMEPREFIX_HUMAN]
                         [--nthreads NTHREADS] [--hit_cutoff HIT_CUTOFF]
                         [-s SUFFIX_SE] [-s1 SUFFIX_PE_1] [-s2 SUFFIX_PE_2]
                         [--fastq]
                         input_path

positional arguments:
  input_path

optional arguments:
  -h, --help            show this help message and exit
  --VIRTUSDir VIRTUSDIR
  --genomeDir_human GENOMEDIR_HUMAN
  --genomeDir_virus GENOMEDIR_VIRUS
  --salmon_index_human SALMON_INDEX_HUMAN
  --salmon_quantdir_human SALMON_QUANTDIR_HUMAN
  --outFileNamePrefix_human OUTFILENAMEPREFIX_HUMAN
  --nthreads NTHREADS
  --hit_cutoff HIT_CUTOFF
  -s SUFFIX_SE, --Suffix_SE SUFFIX_SE
  -s1 SUFFIX_PE_1, --Suffix_PE_1 SUFFIX_PE_1
  -s2 SUFFIX_PE_2, --Suffix_PE_2 SUFFIX_PE_2
  --fastq

[user@cn3144 VIRTUS]$ **cwltool --singularity --tmp-outdir-prefix=/lscratch/$SLURM\_JOB\_ID/ \
 --tmpdir-prefix=/lscratch/$SLURM\_JOB\_ID/ \
 $VIRTUS\_WORKFLOW/VIRTUS.PE.cwl \
 --fastq1 ERR3240275\_1.fastq.gz --fastq2 ERR3240275\_2.fastq.gz \
 --genomeDir\_human $VIRTUS\_INDEX/STAR\_index\_human \
 --genomeDir\_virus $VIRTUS\_INDEX/STAR\_index\_virus \
 --salmon\_index\_human $VIRTUS\_INDEX/salmon\_index\_human \
 --salmon\_quantdir\_human salmon\_human \
 --nthreads 40**
[user@cn3144 ]$ **module load VIRTUS/2.0.1**
[user@cn3144 ]$ **VIRTUS\_wrapper.py -h** 
usage: VIRTUS_wrapper.py [-h] [--VIRTUSDir VIRTUSDIR] --genomeDir_human GENOMEDIR_HUMAN --genomeDir_virus GENOMEDIR_VIRUS
                         [--outFileNamePrefix_human OUTFILENAMEPREFIX_HUMAN] [--nthreads NTHREADS] [-s SUFFIX_SE] [-s1 SUFFIX_PE_1] [-s2 SUFFIX_PE_2]
                         [--fastq] [--figsize FIGSIZE] [--th_cov TH_COV] [--th_rate TH_RATE]
                         input_path

positional arguments:
  input_path

optional arguments:
  -h, --help            show this help message and exit
  --VIRTUSDir VIRTUSDIR
  --genomeDir_human GENOMEDIR_HUMAN
  --genomeDir_virus GENOMEDIR_VIRUS
  --outFileNamePrefix_human OUTFILENAMEPREFIX_HUMAN
  --nthreads NTHREADS
  -s SUFFIX_SE, --Suffix_SE SUFFIX_SE
  -s1 SUFFIX_PE_1, --Suffix_PE_1 SUFFIX_PE_1
  -s2 SUFFIX_PE_2, --Suffix_PE_2 SUFFIX_PE_2
  --fastq
  --figsize FIGSIZE     (default:8,3)
  --th_cov TH_COV       threshold of max viral coverage to plot, test (default:10)
  --th_rate TH_RATE     threshold of max rate virus/human to plot, test (default:0.0001)

[user@cn3144 ]$ **cwltool --singularity \
 --tmp-outdir-prefix=/lscratch/$SLURM\_JOB\_ID/ \
 --tmpdir-prefix=/lscratch/$SLURM\_JOB\_ID/ \
 $VIRTUS\_WORKFLOW/VIRTUS.PE.cwl \
 --fastq1 ERR3240275\_1.fastq.gz \
 --fastq2 ERR3240275\_2.fastq.gz \
 --genomeDir\_human $VIRTUS\_INDEX/STAR\_index\_human \
 --genomeDir\_virus $VIRTUS\_INDEX/STAR\_index\_virus \
 --nthreads 40**
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. VIRTUS.sh). For example:



```

#!/bin/bash
#SBATCH --job-name=S1_VIRTUS
#SBATCH --output=S1_VIRTUS.out
#SBATCH --ntasks=1
#SBATCH --gres=lscratch:500
#SBATCH --cpus-per-task=40
#SBATCH --mem=40Gb
#SBATCH --time=8:00:00
#SBATCH --partition=norm

set -e
module load VIRTUS/2
cp "$VIRTUS_TEST_DATA"/* /data/$USER/VIRTUS
cd /data/$USER/VIRTUS
VIRTUS_wrapper.py input.csv \ 
--genomeDir_human $VIRTUS_INDEX/STAR_index_human \
--genomeDir_virus $VIRTUS_INDEX/STAR_index_virus \
--nthreads 40


```

 Submit the job:

```
sbatch VIRTUS.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.

```


```










