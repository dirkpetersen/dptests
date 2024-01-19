

document.querySelector('title').textContent = 'FreePSI on Biowulf';
FreePSI on Biowulf


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



FreePSI is a new method for genome-wide percent spliced in (PSI) estimation that requires neither a reference transcriptome (hence, transcriptome-free) nor the mapping of RNA-seq reads (hence, alignment-free). The first attribute allows FreePSI to work effectively when a high quality reference transcriptome is unavailable and the second not only helps make FreePSI more efficient, it also eliminates the necessity of dealing with multi-reads.



### References:


* [Zhou, Jianyu, et al. "FreePSI: an alignment-free approach to estimating exon-inclusion ratios without a reference transcriptome." *Nucleic acids research* 46.2 (2018): e11-e11.](https://academic.oup.com/nar/article-abstract/46/2/e11/4607800)


Documentation
* [FreePSI GitHub repo](https://github.com/JY-Zhou/FreePSI)


Important Notes
* Module Name: freepsi (see [the modules page](/apps/modules.html) for more information)
 * Set the number of threads for each instance using the SLURM\_CPUS\_PER\_TASK variable. See examples below.
 * When calling the python scripts, do so without adding the python or python3 directives. That will be added automatically.
* Environment variables set 
	+ FREEPSI\_HOME
	+ FREEPSI\_TESTDATA* Example files in /usr/local/apps/freepsi/TESTDATA (pointed to by the FREEPSI\_TESTDATA environment variable).



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. In this example, we run a script that contains FreePSI commands. Reviewing the script provides individual command syntax.
  
Sample session (user input in **bold**):



```

[user@biowulf ~]$ **sinteractive -c6 --mem=12g --gres=lscratch:10**
salloc: Pending job allocation 33604907
salloc: job 33604907 queued and waiting for resources
salloc: job 33604907 has been allocated resources
salloc: Granted job allocation 33604907
salloc: Waiting for resource configuration
salloc: Nodes cn0849 are ready for job
srun: error: x11: no local DISPLAY defined, skipping
error: unable to open file /tmp/slurm-spank-x11.33604907.0
slurmstepd: error: x11: unable to read DISPLAY value

[user@cn0849 ~]$ **module load freepsi**
[+] Loading freepsi  0.3  on cn0849
[+] Loading singularity  3.8.5-1  on cn0849
[+] Loading jellyfish  2.3.0

[user@cn0849 ~]$ **cd /lscratch/$SLURM\_JOB\_ID**

[user@cn0849 33604907]$ **cp -r $FREEPSI\_TESTDATA .**

[user@cn0849 33604907]$ **cd TESTDATA/**

[user@cn0849 TESTDATA]$ **cat run\_FreePSI.sh**
#!/bin/bash

set -e

# Provide the directory containing Jellyfish (usually named as 'bin')
Jellyfish="$(dirname $(which jellyfish))"
# Provide the directory containing FreePSI
# E.g.
#FreePSI=../bin
FreePSI="${FREEPSI_HOME}/bin"

if [ -z ${Jellyfish} ]; then
    echo "Error: Please modify me to provide the directory containing Jellyfish"
    exit;
fi

if [ -z ${FreePSI} ]; then
    echo "Error: Please modify me to provide the directory containing FreePSI"
    exit;
fi

GENOME_DIR=./genome
BND_FILE=./annotation/hg38_refGene_exonBoundary_chr21.bed
READS=RNA-seq
K=27
THREAD=$SLURM_CPUS_PER_TASK

set -x
# Count k-mers in RNA-seq reads using jellyfish
${Jellyfish}/jellyfish count -m ${K} -s 100M -t ${THREAD} -Q 5 ${READS}/reads_final.1.fastq -o ${READS}/reads.1.jf
${Jellyfish}/jellyfish dump ${READS}/reads.1.jf -o ${READS}/reads.1.fa
${Jellyfish}/jellyfish count -m ${K} -s 100M -t ${THREAD} -Q 5 ${READS}/reads_final.2.fastq -o ${READS}/reads.2.jf
${Jellyfish}/jellyfish dump ${READS}/reads.2.jf -o ${READS}/reads.2.fa

# Produce raw estimates of PSI values using FreePSI
#Build
${FreePSI}/freePSI build\
    -k $K -p ${THREAD} \
    -g ${GENOME_DIR} \
    -a ${BND_FILE} \
    -1 ${READS}/reads.1.fa \
    -2 ${READS}/reads.2.fa \
    -o ./hashtable.json

#Quant
${FreePSI}/freePSI quant\
    -k $K -p ${THREAD} \
    -i ./hashtable.json \
    -o .

# Post-process the raw estimates of PSI values
# python3 ${FreePSI}/postProc.py \
${FreePSI}/postProc.py \
    ./psi_freePSI_raw.json \
    ./psi_freePSI.json

# Summarize the PSI values into a readable file
# python3 ${FreePSI}/summary.py \
${FreePSI}/summary.py \
    ./annotation/hg38_refGene_exonBoundary_chr21.bed \
    ./psi_freePSI.json \
    ./psi_freePSI.summary

[user@cn0849 TESTDATA]$ **./run\_FreePSI.sh**
+ /usr/local/apps/jellyfish/2.3.0/bin/jellyfish count -m 27 -s 100M -t 6 -Q 5 RNA-seq/reads_final.1.fastq -o RNA-seq/reads.1.jf
+ /usr/local/apps/jellyfish/2.3.0/bin/jellyfish dump RNA-seq/reads.1.jf -o RNA-seq/reads.1.fa
+ /usr/local/apps/jellyfish/2.3.0/bin/jellyfish count -m 27 -s 100M -t 6 -Q 5 RNA-seq/reads_final.2.fastq -o RNA-seq/reads.2.jf
+ /usr/local/apps/jellyfish/2.3.0/bin/jellyfish dump RNA-seq/reads.2.jf -o RNA-seq/reads.2.fa
+ /usr/local/apps/freepsi/0.3/bin/freePSI build -k 27 -p 6 -g ./genome -a ./annotation/hg38_refGene_exonBoundary_chr21.bed -1 RNA-seq/reads.1.fa -2 RNA-seq/reads.2.fa -o ./hashtable.json

### Start to build theoretical and real kmer profile ...
[...snip]
Elasped time 2s.

### Start to refine solution and compute PSI ...
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> OK!
### Finish refining solution and computing PSI ...
### Finished!
CPU Time elapsed: 7s.
Natural Time elapsed: 3s.
+ /usr/local/apps/freepsi/0.3/bin/postProc.py ./psi_freePSI_raw.json ./psi_freePSI.json
+ /usr/local/apps/freepsi/0.3/bin/summary.py ./annotation/hg38_refGene_exonBoundary_chr21.bed ./psi_freePSI.json ./psi_freePSI.summary

[user@cn0849 TESTDATA]$ **exit**
exit
salloc: Relinquishing job allocation 33604907

[user@biowulf ~]$ 



```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. freepsi.sh). Use the script in TESTDATA (shown above) as an example:


Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] freepsi.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. freepsi.swarm). For the sake of simplicity we place a single command in each line, but you can place multiple commands on a single line seperated by semicolons or you can create scripts and execute a different script or the same script with different inputs on each line instead:



```

freePSI quant -k $K -p ${SLURM_CPUS_PER_TASK} -i /path/to/hashtable1.json -o /our/dir1
freePSI quant -k $K -p ${SLURM_CPUS_PER_TASK} -i /path/to/hashtable2.json -o /our/dir2
freePSI quant -k $K -p ${SLURM_CPUS_PER_TASK} -i /path/to/hashtable3.json -o /our/dir3
freePSI quant -k $K -p ${SLURM_CPUS_PER_TASK} -i /path/to/hashtable4.json -o /our/dir4

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f freepsi.swarm [-g #] [-t #] --module freepsi
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module freepsi Loads the freepsi module for each subjob in the swarm 
 | |
 | |
 | |








