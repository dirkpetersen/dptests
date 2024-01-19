

document.querySelector('title').textContent = 'adVNTR on Biowulf';
adVNTR on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
 |


 adVNTR is a tool for genotyping Variable Number Tandem Repeats (VNTR) from sequence data. It works with both NGS short reads (Illumina HiSeq) and SMRT reads (PacBio) and finds diploid repeating counts for VNTRs and identifies possible mutations in the VNTR sequences.



### References:


* Bakhtiari, M., Shleizer-Burko, S., Gymrek, M., Bansal, V. and Bafna, V.
*Targeted genotyping of variable number tandem repeats with adVNTR.* Genome Research, 28(11), pp.1709-1719.
[Journal](https://genome.cshlp.org/content/28/11/1709/)


Documentation
* advntr on [GitHub](https://github.com/mehrdadbakhtiari/adVNTR)


Important Notes
* Module Name: advntr (see [the modules page](/apps/modules.html)
for more information)
* Example files in `$ADVNTR_TEST_DATA`
* The warning messages at the top regarding cudart library when running `advntr` can be safely ignored.



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=5g --gres=lscratch:20**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144]$ **cd /lscratch/$SLURM\_JOB\_ID**
[user@cn3144]$ **module load advntr**
[user@cn3144]$ **cp ${ADVNTR\_TEST\_DATA} .**
[user@cn3144]$ **cd TEST\_DATA**
[user@cn3144]$ **mkdir log\_dir**
[user@cn3144]$ **advntr genotype --vntr\_id 301645 --alignment\_file CSTB\_2\_5\_testdata.bam --working\_directory log\_dir/**
    2021-04-30 18:46:02.676445: I tensorflow/stream_executor/platform/default/dso_loader.cc:49] Successfully opened dynamic library libcudart.so.11.0
    301645
    2/2

```



Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch script file (e.g. advntr.sh). For example:



```

#!/bin/bash
cd /lscratch/$SLURM_JOB_ID
module load advntr
cp $ADVNTR_TEST_DATA .
cd TEST_DATA
advntr genotype --vntr_id 301645 --alignment_file CSTB_2_5_testdata.bam --working_directory log_dir/
....
....

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --mem=10g --gres=lscratch:20 advntr.sh
```







