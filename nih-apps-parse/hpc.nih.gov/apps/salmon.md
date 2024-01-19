

document.querySelector('title').textContent = 'salmon on Biowulf';
salmon on Biowulf


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


[![Salmon Logo](/images/SalmonLogo.png)](http://dasher.wustl.edu)

Salmon is used to estimate expression at the transcript level from RNA-Seq data. It uses
*quasi-mapping* instead of full alignment which reduces computation time and
attempts to account for common biases in RNA-Seq data with realistic models.



### References:


* R. Patro, G. Duggal, M. Love, R. Irizarry, and C. Kingsford.
 *Salmon provides fast and bias-aware quantification of transcript expression*.
 Nature Methods 2017, 14:417-419.
 [PubMed](https://www.ncbi.nlm.nih.gov/pubmed/28263959) | 
 [PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5600148/) | 
 [Journal](https://www.nature.com/articles/nmeth.4197)


Documentation
* Salmon [Home page](https://combine-lab.github.io/salmon/about/)
* Salmon [GitHub repo](https://github.com/COMBINE-lab/salmon)


Important Notes
* Module Name: salmon (see [the modules page](/apps/modules.html) for more information)
* Salmon is a multi-threaded application. Please make sure to match the number of threads to the
 number of allocated CPUs.
* Example files can be found in `$SALMON_TEST_DATA`
* When running large numbers of concurrent jobs with a large transcriptome index,
 please consider copying the transcriptome index to lscratch and running several (5-10) consecutive
 salmon analyses with the same local index to reduce the load on the file systems.



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.
In the sample session below, we will build a transcriptome index for *S. cerevisiae* and then
use that index to quantify expression in some yeast RNA-Seq data:



```

[user@biowulf]$ **sinteractive --cpus-per-task=6 --gres=lscratch:10**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load salmon**
[user@cn3144 ~]$ **cd /lscratch/$SLURM\_JOB\_ID**
[user@cn3144 ~]$ **cp -rL $SALMON\_TEST\_DATA/\* .**
[user@cn3144 ~]$ **ls -lh**
total 9.0M
drwxr-xr-x 2 user group 4.0K Feb 23 07:34 ERP004763
-rw-r--r-- 1 user group 9.0M Feb 23 07:34 R64-1-1.cdna_nc.fa

[user@cn3144 ~]$ # create the transcriptome index from the fasta of transcripts
[user@cn3144 ~]$ **salmon index -p $SLURM\_CPUS\_PER\_TASK -t R64-1-1.cdna\_nc.fa -i R64**
[2020-04-22 11:37:08.511] [jLog] [warning] The salmon index is being built without any decoy sequences.  It is recommended that decoy sequence (either computed auxiliary decoy sequence or the genome of the organism) be provided during indexing. Further details can be found at https://salmon.readthedocs.io/en/latest/salmon.html#preparing-transcriptome-indices-mapping-based-mode.
[2020-04-22 11:37:08.511] [jLog] [info] building index
[...snip...]
[2020-04-22 11:37:13.870] [jLog] [info] done building index
[user@cn3144 ~]$ # quantify the 6 samples in ERP004763 using automatically determined library
                 # type (e.g. stranded vs unstranded; -l A)
[user@cn3144 ~]$ **for fq in ERP004763/\*.fastq.gz; do
 echo "Processing $fq"
 salmon quant -i R64 -l A -r $fq -p $SLURM\_CPUS\_PER\_TASK \
 -o $(basename $fq .fastq.gz)\_quant
 done** 
[...snip...]
[user@cn3144 ~]$ **ls -lh ERR458495\_quant**
total 296K
drwxr-xr-x 2 user group 4.0K Feb 23 07:43 aux_info
-rw-r--r-- 1 user group  204 Feb 23 07:43 cmd_info.json
-rw-r--r-- 1 user group  490 Feb 23 07:43 lib_format_counts.json
drwxr-xr-x 2 user group 4.0K Feb 23 07:43 libParams
drwxr-xr-x 2 user group 4.0K Feb 23 07:43 logs
-rw-r--r-- 1 user group 272K Feb 23 07:43 quant.sf

[user@cn3144 ~]$ **head -5 ERR458495\_quant/quant.sf**
Name    Length  EffectiveLength TPM     NumReads
YHR055C 186     8.490   2588.845789     88.963768
YPR161C 1974    1725.000        6.445037        45.000000
YOL138C 4026    3777.000        3.008935        46.000000
YDR395W 2835    2586.000        7.834068        82.000000

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. salmon.sh) similar to the following example:



```

#!/bin/bash

module load salmon/0.9.1 || exit 1

fastq=${1:-none}
[[ "$fastq" != "none" ]] || exit 1
[[ -f "$fastq" ]] || exit 1

salmon quant -i R64 -l A -r $fastq -p $SLURM_CPUS_PER_TASK \
    -o $(basename $fastq .fastq.gz)_quant

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=6 --mem=4g salmon.sh ERP004763/ERR458495.fastq.gz
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. salmon.swarm). For example:



```

salmon quant -i R64 -l A -r sample1.fastq.gz -p $SLURM_CPUS_PER_TASK -o sample1_quant
salmon quant -i R64 -l A -r sample2.fastq.gz -p $SLURM_CPUS_PER_TASK -o sample2_quant
salmon quant -i R64 -l A -r sample3.fastq.gz -p $SLURM_CPUS_PER_TASK -o sample3_quant

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f salmon.swarm -g 4 -t 6 --module salmon/0.9.1
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g #  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t #  Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module salmon  Loads the salmon module for each subjob in the swarm
 | |
 | |
 | |








