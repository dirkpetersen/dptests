

document.querySelector('title').textContent = 'mantis on Biowulf';
mantis on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
 |


MANTIS (Microsatellite Analysis for Normal-Tumor InStability) is a program developed for detecting microsatellite instability from paired-end BAM files. To perform analysis, the program needs a tumor BAM and a matched normal BAM file (produced using the same pipeline) to determine the instability score between the two samples within the pair. Longer reads (ideally, 100 bp or longer) are recommended, as shorter reads are unlikely to entirely cover the microsatellite loci, and will be discarded after failing the quality control filters.



### References:


* Kautto, E. A., Bonneville, R., Miya, J., Yu, L., Krook, M. A., Reeser, J. W., & Roychowdhury, S
*Performance evaluation for rapid detection of pan-cancer microsatellite instability with MANTIS.*Oncotarget, 8(5), 7452–7463.
<http://doi.org/10.18632/oncotarget.13918>


Documentation
* mantis on [GitHub](https://github.com/OSU-SRLab/MANTIS)


Important Notes
* Module Name: mantis (see [the modules page](/apps/modules.html)
for more information)
* mantis is multithreaded. Please match the number of threads to the number of allocated CPUs
* Example files in `${MANTIS_TEST_DATA}`



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=5g --cpus-per-task=4**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144]$ **cd /data/$USER**
[user@cn3144]$ **cp ${MANTIS\_TEST\_DATA} .**
[user@cn3144]$ **cd TEST\_DATA**
[user@cn3144]$ **HG19=/fdb/igenomes/Homo\_sapiens/UCSC/hg19/hg19.fa**
[user@cn3144]$ **module load mantis**
[user@cn3144]$ **mantis.py --threads $SLURM\_CPUS\_PER\_TASK --bedfile target\_loci.bed --genome $HG19 -n normal.bam -t tumor.bam -mrq 20.0 -mlq 25.0 -mlc 20 -mrr 1 -o mantis.out**
    

```



Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch script file (e.g. mantis.sh). For example:



```

#!/bin/bash
module load mantis
cd /data/$USER
cp -r ${MANTIS_TEST_DATA} .
cd TEST_DATA
mantis.py --threads $SLURM_CPUS_PER_TASK --bedfile target_loci.bed --genome $HG19 -n normal.bam -t tumor.bam -mrq 20.0 -mlq 25.0 -mlc 20 -mrr 1 -o mantis.out

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --mem=5g --cpus-per-task=4 mantis.sh
```







