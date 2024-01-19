

document.querySelector('title').textContent = 'kneaddata on Biowulf';
kneaddata on Biowulf


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



From the KneadData user manual:




>  
>  KneadData is a tool designed to perform quality control on
>  metagenomic and metatranscriptomic sequencing data, especially data from
>  microbiome experiments. In these experiments, samples are typically taken
>  from a host in hopes of learning something about the microbial community on
>  the host. However, sequencing data from such experiments will often contain
>  a high ratio of host to bacterial reads. This tool aims to perform
>  principled in silico separation of bacterial reads from these "contaminant"
>  reads, be they from the host, from bacterial 16S sequences, or other
>  user-defined sources. Additionally, KneadData can be used for other
>  filtering tasks. For example, if one is trying to clean data derived from a
>  human sequencing experiment, KneadData can be used to separate the human
>  and the non-human reads. 
> 


Documentation
* [Manual](https://huttenhower.sph.harvard.edu/kneaddata/)


Important Notes
* Module Name: kneaddata (see [the modules page](/apps/modules.html) for more information)
* kneaddata is a multithreaded/multiprocess application. Make sure to match the number 
 of cpus requested with the number of threads.
* Example files in `$KNEADDATA_TEST_DATA`
* Reference data in /fdb/kneaddata/



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:


Allocate an interactive session with [sinteractive](/docs/userguide.html#int)
and use as shown below. In this case we will use test data that is 
an artificial mixture of 1M human exome reads and 1M environmental metagenomic 
reads. The 50% human reads is treated as an artificial contamination and removed:



```

[user@biowulf]$ **sinteractive --mem=12g --cpus-per-task=8 --gres=lscratch:10**
salloc.exe: Pending job allocation 33247354
salloc.exe: job 33247354 queued and waiting for resources
salloc.exe: job 33247354 has been allocated resources
salloc.exe: Granted job allocation 33247354
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job
srun: error: x11: no local DISPLAY defined, skipping

[user@nc3144]$ **module load kneaddata**
[user@nc3144]$ **cd /lscratch/${SLURM\_JOB\_ID}**
[user@nc3144]$ # copy test data
[user@nc3144]$ **cp $KNEADDATA\_TEST\_DATA/\* .**
[user@nc3144]$ # the test data is 50% human, 50% environmental metagenomic data
[user@nc3144]$ # and the read names are labelled accordingly
[user@nc3144]$ **zcat test\_R1.fastq.gz \
 | awk '/@meta/ {m++} /@human/ {h++} END {printf("human: %i\nmeta: %i\n", h, m)}'**
human: 1000000
meta:  1000000
[user@nc3144]$ **# run kneaddata on the paired end test data**
[user@nc3144]$ **kneaddata -i test\_R1.fastq.gz -i test\_R2.fastq.gz \
 --reference-db $KNEADDATA\_DB/Homo\_sapiens\_Bowtie2\_v0.1/Homo\_sapiens \
 --output-prefix test --output test\_out \
 -p 2 -t 4 --run-fastqc-end** 
[...snip...]
Final output files created: 
/lscratch/33247583/test_out/test_paired_1.fastq
/lscratch/33247583/test_out/test_paired_2.fastq
/lscratch/33247583/test_out/test_unmatched_1.fastq
/lscratch/33247583/test_out/test_unmatched_2.fastq

[user@nc3144]$ # check the composition of human/metagenome reads in the cleand data
[user@nc3144]$ **cat test\_out/test\_paired\_1.fastq \
 | awk '/@meta/ {m++} /@human/ {h++} END {printf("human: %i\nmeta: %i\n", h, m)}'**
human: 18641
meta:  951886
[user@nc3144]$ **exit**
[user@biowulf]$

```

So the 50% artificial contamination with human reads was reduced to 2%.



Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. kneaddata.sh) similar to the following example:



```

#! /bin/bash

module load kneaddata/0.7.0 || exit 1
cd /lscratch/${SLURM_JOB_ID} || exit 1

if [[ ! -e test_R1.fastq.gz ]]; then
    cp $KNEADDATA_TEST_DATA/* .
fi
rm -rf test_out

kneaddata -i test_R1.fastq.gz -i test_R2.fastq.gz \
  --reference-db $KNEADDATA_DB/Homo_sapiens_Bowtie2_v0.1/Homo_sapiens \
  --output-prefix test --output test_out \
  -p 2 -t 4 --run-fastqc-end

cp -r test_out /data/$USER/important_project

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=8 --mem=12g --gres=lscratch:10 kneaddata.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. kneaddata.swarm). For example:



```

kneaddata -i sample1_R1.fastq.gz -i sample1_R2.fastq.gz \
  --reference-db $KNEADDATA_DB/Homo_sapiens_Bowtie2_v0.1/Homo_sapiens \
  --output-prefix sample1 --output sample1_out -p 2 -t 4 --run-fastqc-end
kneaddata -i sample2_R1.fastq.gz -i sample2_R2.fastq.gz \
  --reference-db $KNEADDATA_DB/Homo_sapiens_Bowtie2_v0.1/Homo_sapiens \
  --output-prefix sample2 --output sample2_out -p 2 -t 4 --run-fastqc-end
kneaddata -i sample3_R1.fastq.gz -i sample2_R2.fastq.gz \
  --reference-db $KNEADDATA_DB/Homo_sapiens_Bowtie2_v0.1/Homo_sapiens \
  --output-prefix sample3 --output sample3_out -p 2 -t 4 --run-fastqc-end

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f kneaddata.swarm -g 12 -t 8 --module kneaddata/0.7.0
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g #  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t #  Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module kneaddata  Loads the kneaddata module for each subjob in the swarm 
 | |
 | |
 | |








