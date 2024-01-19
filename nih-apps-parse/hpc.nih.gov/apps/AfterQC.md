

document.querySelector('title').textContent = 'Afterqc on Biowulf';
Afterqc on Biowulf


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



Automatic Filtering, Trimming, Error Removing and Quality Control for fastq data. AfterQC can simply go through all fastq files in a folder and then output three folders: good, bad and QC folders, which contains good reads, bad reads and the QC results of each fastq file/pair.



### References:


* Shifu Chen, Tanxiao Huang, Yanqing Zhou, Yue Han, Mingyan Xu and Jia Gu. AfterQC: automatic filtering, trimming, error removing and quality control for fastq data. BMC Bioinformatics 2017 18(Suppl 3):80 https://doi.org/10.1186/s12859-017-1469-3


Documentation
* [Afterqc Main Site](https://github.com/OpenGene/AfterQC)


Important Notes
* Module Name: afterqc (see [the modules page](/apps/modules.html) for more information)
* Single threaded app
* environment variables set 
	+ afterqc\_HOME=/usr/local/apps/afterqc/$version



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session (user input in **bold**:



```

[user@biowulf ~]$ **sinteractive -c4 --mem=8g --gres=lscratch:10**
salloc.exe: Pending job allocation 11086981
salloc.exe: job 11086981 queued and waiting for resources
salloc.exe: job 11086981 has been allocated resources
salloc.exe: Granted job allocation 11086981
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn0882 are ready for job
srun: error: x11: no local DISPLAY defined, skipping
error: unable to open file /tmp/slurm-spank-x11.11086981.0
slurmstepd: error: x11: unable to read DISPLAY value

[user@cn0882 ~]$ **cd /lscratch/$SLURM\_JOB\_ID**

[user@cn0882 11086981]$ **module load afterqc**
[+] Loading afterqc  0.9.7  on cn0882
[+] Loading python 3.7  ...

[user@cn0882 11086981]$ **python2 $afterqc\_HOME/after.py \
 -1 /fdb/app\_testdata/fastq/H\_sapiens/hg100\_1m\_pe1.fq.gz \
 -2 /fdb/app\_testdata/fastq/H\_sapiens/hg100\_1m\_pe2.fq.gz**

/fdb/app_testdata/fastq/H_sapiens/hg100_1m_pe1.fq.gz options:
{'qc_only': False, 'version': '0.9.6', 'seq_len_req': 35, 'index1_file': None, 'trim_tail': 0, 'report_output_folder': None, 'trim_pair_same': True, 
'no_correction': False, 'debubble_dir': 'debubble', 'barcode_flag': 'barcode', 'read2_file': '/fdb/app_testdata/fastq/H_sapiens/hg100_1m_pe2.fq.gz', 
'barcode_length': 12, 'trim_tail2': 0, 'unqualified_base_limit': 60, 'allow_mismatch_in_poly': 2, 'read2_flag': 'R2', 'store_overlap': False, 
'debubble': False, 'read1_flag': 'R1', 'index2_flag': 'I2', 'draw': True, 'index1_flag': 'I1', 'mask_mismatch': False, 'barcode': False, 
'gzip': False, 'overlap_output_folder': None, 'barcode_verify': 'CAGTA', 'compression': 2, 'index2_file': None, 'qualified_quality_phred': 15, 
'trim_front': 1, 'good_output_folder': 'good', 'poly_size_limit': 35, 'n_base_limit': 5, 'qc_sample': 200000, 'trim_front2': 1, 'no_overlap': False, 
'input_dir': None, 'read1_file': '/fdb/app_testdata/fastq/H_sapiens/hg100_1m_pe1.fq.gz', 'qc_kmer': 8, 'bad_output_folder': None}

Time used: 544.626863003

[user@cn0882 11086981]$ **ls -lh**
total 12K
drwxr-xr-x 2 user user 4.0K Mar 22 17:47 bad
drwxr-xr-x 2 user user 4.0K Mar 22 17:47 good
drwxr-xr-x 2 user user 4.0K Mar 22 17:54 QC

[user@cn0882 11086981]$ **exit**
exit
salloc.exe: Relinquishing job allocation 11086981
salloc.exe: Job allocation 11086981 has been revoked.

[user@biowulf ~]$ 

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. afterqc.sh). For example:



```

#!/bin/bash
set -e
module load afterqc
python2 $afterqc_HOME/after.py -1 /fdb/app_testdata/fastq/H_sapiens/hg100_1m_pe1.fq.gz -2 /fdb/app_testdata/fastq/H_sapiens/hg100_1m_pe2.fq.gz 

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch afterqc.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. afterqc.swarm). For example:



```

python2 $afterqc_HOME/after.py -1 sample1_1.fq.gz -2 sample1_2.fq.gz 
python2 $afterqc_HOME/after.py -1 sample2_1.fq.gz -2 sample2_2.fq.gz 
python2 $afterqc_HOME/after.py -1 sample3_1.fq.gz -2 sample3_2.fq.gz 

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f afterqc.swarm --module afterqc
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module afterqc Loads the afterqc module for each subjob in the swarm 
 | |
 | |
 | |








