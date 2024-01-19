

document.querySelector('title').textContent = 'SpliceAI: predicting splicing from primary sequence with deep learning';
**SpliceAI: predicting splicing from primary sequence with deep learning**


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


  


SpliceAI is a deep neural network that accurately predicts
splice junctions from an arbitrary pre-mRNA
transcript sequence, enabling precise prediction of
noncoding genetic variants that cause cryptic
splicing.



  

### References:


* Kishore Jaganathan, Sofia Kyriazopoulou Panagiotopoulou, Jeremy F. McRae, et al.   

*Predicting Splicing from Primary Sequence with Deep Learning* ,
Cell **176**, 535–548, January 24, 2019


Documentation
* [SpliceAI GitHub page](https://github.com/Illumina/SpliceAI)


Important Notes
* Module Name: SpliceAI (see [the modules page](/apps/modules.html) for more information)
* Singlethreaded
* Unusual environment variables set
	+ **SPLICEAI\_HOME**   installation directory
	+ **SPLICEAI\_BIN**         executable directory
	+ **SPLICEAI\_DATA**     sample data directory
	+ **SPLICEAI\_SRC**       source code directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=8g --gres=gpu:v100:1**
[user@cn4466 ~]$ **module load SpliceAI**
+] Loading python 3.6  ...
[+] Loading cuDNN 7.0  libraries...
[+] Loading CUDA Toolkit  9.0.176  ...
[user@cn4466 ~]$ **spliceai -h**
...
usage: spliceai [-h] [-I [input]] [-O [output]] -R reference -A annotation
                [-D [distance]] [-M [mask]]

Version: 1.3

optional arguments:
  -h, --help     show this help message and exit
  -I [input]     path to the input VCF file, defaults to standard in
  -O [output]    path to the output VCF file, defaults to standard out
  -R reference   path to the reference genome fasta file
  -A annotation  "grch37" (GENCODE V24lift37 canonical annotation file in
                 package), "grch38" (GENCODE V24 canonical annotation file in
                 package), or path to a similar custom gene annotation file
  -D [distance]  maximum distance between the variant and gained/lost splice
                 site, defaults to 50
  -M [mask]      mask scores representing annotated acceptor/donor gain and
                 unannotated acceptor/donor loss, defaults to 0

```

Download sample data:

```

[user@cn4466 ~]$ **cp $SPLICEAI\_DATA/\* .** 

```

Specify a reference sequence:

```

[user@cn4466 ~]$ **ln -s /fdb/igenomes/Homo\_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa hg38.fa**

```

Run the spliceai executable on the sample data:

```
 
[user@cn4466 ~]$ **spliceai -I input.vcf -O output.vcf -R hg38.fa -A grch37 &** 
...
[user@cn4466 ~]$ **nvidia-smi**
...
+-----------------------------------------------------------------------------+
| NVIDIA-SMI 470.82.01    Driver Version: 470.82.01    CUDA Version: 11.4     |
|-------------------------------+----------------------+----------------------+
| GPU  Name        Persistence-M| Bus-Id        Disp.A | Volatile Uncorr. ECC |
| Fan  Temp  Perf  Pwr:Usage/Cap|         Memory-Usage | GPU-Util  Compute M. |
|                               |                      |               MIG M. |
|===============================+======================+======================|
|   0  Tesla K80           On   | 00000000:84:00.0 Off |                  Off |
| N/A   35C    P0    71W / 149W |  11621MiB / 12206MiB |      0%      Default |
|                               |                      |                  N/A |
+-------------------------------+----------------------+----------------------+

+-----------------------------------------------------------------------------+
| Processes:                                                                  |
|  GPU   GI   CI        PID   Type   Process name                  GPU Memory |
|        ID   ID                                                   Usage      |
|=============================================================================|
|    0   N/A  N/A     45210      C   /usr/bin/python                 11616MiB |
+-----------------------------------------------------------------------------+
...
WARNING:root:Skipping record (ref issue): 2     152389953       .       T       A,C,G   .       .       .

WARNING:root:Skipping record (ref issue): 2     179415988       .       C       CA      .       .       .

WARNING:root:Skipping record (ref issue): 2     179446218       .       ATACT   A       .       .       .

WARNING:root:Skipping record (ref issue): 2     179446218       .       ATACT   AT,ATA  .       .       .

WARNING:root:Skipping record (ref issue): 2     179642185       .       G       A       .       .       .

WARNING:root:Skipping record (ref issue): 19    38958362        .       C       T       .       .       .

WARNING:root:Skipping record (ref issue): 21    47406854        .       CCA     C       .       .       .

WARNING:root:Skipping record (ref issue): 21    47406856        .       A       AT      .       .       .

WARNING:root:Skipping record (ref issue): X     129274636       .       A       C,G,T   .       .       .

```

An output file output.vcf will be produced:

```

[user@cn4466 ~]$ **cat output.vcf**
...
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
1       25000   .       A       C,G,T   .       .       .
2       152389953       .       T       A,C,G   .       .       .
2       179415988       .       C       CA      .       .       .
2       179446218       .       ATACT   A       .       .       .
2       179446218       .       ATACT   AT,ATA  .       .       .
2       179642185       .       G       A       .       .       .
19      38958362        .       C       T       .       .       .
21      47406854        .       CCA     C       .       .       .
21      47406856        .       A       AT      .       .       .
X       129274636       .       A       C,G,T   .       .       .

```

Exit the application:

```

[user@cn4466 ~]$ **exit**
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. spliceai.sh). For example:



```

#!/bin/bash
module load SpliceAI      
cp $SPLICEAI_DATA/* .  
spliceai -I input.vcf -O output.vcf -R /fdb/igenomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa  -A grch37

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] spliceai.sh
```







