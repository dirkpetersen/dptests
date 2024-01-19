

document.querySelector('title').textContent = 'Metal on HPC';
Metal on HPC


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


 The METAL software is designed to facilitate meta-analysis of large datasets 
 (such as several whole genome scans) in a convenient, rapid and memory efficient 
 manner. 


### 


Documentation
* <https://github.com/statgen/METAL>



Important Notes
* Module Name: metal (see [the modules 
 page](/apps/modules.html) for more information)
* Example files in /usr/local/apps/metal/examples/





Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load metal**
[user@cn3144 ~]$ **cd /data/$USER/metal/examples/GlucoseExample**
[user@cn3144 ~]$ **metal << EOT**   
MARKER SNP  
WEIGHT N  
ALLELE EFFECT_ALLELE NON_EFFECT_ALLELE  
FREQ EFFECT_ALLELE_FREQ  
EFFECT BETA  
STDERR SE  
PVAL P_VAL
PROCESS DGI\_three\_regions.txt

MARKER SNP  
 ALLELE EFFECT\_ALLELE NON\_EFFECT\_ALLELE  
 FREQ FREQ\_EFFECT  
 WEIGHT N  
 EFFECT BETA  
 STDERR SE  
 PVAL PVALUE

PROCESS MAGIC\_FUSION\_Results.txt

MARKER SNP  
 DEFAULT 4108  
 ALLELE AL1 AL2  
 FREQ FREQ1  
 EFFECT EFFECT  
 STDERR SE  
 PVAL PVALUE

PROCESS magic\_SARDINIA.tbl

ANALYZE

Press Ctrl-D to submit these commands. Two output files will be created.
[...etc...]

[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$



```






