

document.querySelector('title').textContent = 'msisensor-pro on Biowulf';
msisensor-pro on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
 |


MSIsensor-pro is an updated version of msisensor. MSIsensor-pro evaluates Microsatellite Instability (MSI) for cancer patients with next generation sequencing data. It accepts the whole genome sequencing, whole exome sequencing and target region (panel) sequencing data as input. MSIsensor-pro introduces a multinomial distribution model to quantify polymerase slippages for each tumor sample and a discriminative sites selection method to enable MSI detection without matched normal samples. For samples of various sequencing depths and tumor purities, MSIsensor-pro significantly outperformed the current leading methods in terms of both accuracy and computational cost. If you want to know more detail about MSIsensor-pro, please see the [MSIsensor-pro Schematics and Internals MSIsensor-pro](https://github.com/xjtu-omics/msisensor-pro/wiki/MSIsensor-pro-Schematics-and-Internals) page.



### References:


* Peng Jia, Xiaofei Yang, Li Guo, Bowen Liu, Jiadong Lin, Hao Liang, et al.
*MSIsensor-pro: fast, accurate, and matched-normal-sample-free detection of microsatellite instability.*Genomics Proteomics Bioinformatics 2020,18(1).
[PDF](https://www.sciencedirect.com/science/article/pii/S1672022920300218)


Documentation
* msisensor-pro on [GitHub](https://github.com/xjtu-omics/msisensor-pro)


Important Notes
* Module Name: msisensor-pro (see [the modules page](/apps/modules.html)
for more information)
* Example files in `${MSISENSORPRO_TEST_DATA}`



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=5g**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

module load msisensor-pro
[user@cn3144]$ **cp -rp ${MSISENSORPRO\_TEST\_DATA} /data/$USER**
[user@cn3144]$ **cd /data/$USER/demo/scripts**
[user@cn3144]$ **./1\_test\_scan.sh**
[user@cn3144]$ **./2\_test\_msi.sh**
[user@cn3144]$ **./3\_test\_baseline.sh**
[user@cn3144]$ **./4\_test\_pro.sh**
    

```



Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch script file (e.g. msisensor-pro.sh). For example:



```

#!/bin/bash
module load msisensor-pro
cp -rp ${MSISENSORPRO_TEST_DATA} /data/$USER
cd /data/$USER/demo/scripts
./1_test_scan.sh
./2_test_msi.sh
./3_test_baseline.sh
./4_test_pro.sh

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --mem=5g msisensor-pro.sh
```







