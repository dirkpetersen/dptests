

document.querySelector('title').textContent = 'multiqc on Biowulf';
multiqc on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
 |


 MultiQC searches a directory tree for summary output or logs of many well
known bioinformatics tools and compiles them into a single html report across
all detected samples and tools. 


The tools MultiQC knows about fall into three groups:



**Pre-alignment tools**
FastQC, Skewer, Cutadapt, ...
**Alignment tools**
Bowtie1, Bowtie2, STAR, ...
**Post-alignment tools**
Samtools, Samblaster, Picard, featureCounts, ...

### References:


* P. Ewels, M. Magnusson, S. Lundin and M. Käller. *MultiQC: Summarize analysis results for multiple tools and samples in a single report*. Bioinformatics 2016.
 [PubMed](https://www.ncbi.nlm.nih.gov/pubmed/27312411) | 
 [PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5039924/) | 
 [Journal](https://academic.oup.com/bioinformatics/article/32/19/3047/2196507)


Documentation
* [Manual](http://multiqc.info/)
* [GitHub](https://www.github.com/ewels/MultiQC)


Important Notes
* Module Name: multiqc (see [the modules page](/apps/modules.html) for more information)
* Example files can be found in `$MULTIQC_TEST_DATA`



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

[user@cn3144 ~]$ **module load multiqc**
[+] Loading multiqc  1.9
[user@cn3144 ~]$ **cp -r $MULTIQC\_TEST\_DATA/wgs .**
[user@cn3144 ~]$ **multiqc --title 'WGS\_test' wgs**
[INFO   ]         multiqc : This is MultiQC v1.9
[INFO   ]         multiqc : Template    : default
[INFO   ]         multiqc : Report title: WGS_test
[INFO   ]         multiqc : Searching   : wgs
Searching 1346 files..  [####################################]  100%
[INFO   ]        qualimap : Found 12 BamQC reports
[INFO   ]          snpeff : Found 6 reports
[INFO   ]     varianteval : Found 6 VariantEval reports
[INFO   ]          picard : Found 6 MarkDuplicates reports
[INFO   ]    fastq_screen : Found 12 reports
[INFO   ]          fastqc : Found 12 reports
[INFO   ]         multiqc : Compressing plot data
[INFO   ]         multiqc : Report      : WGS_test_multiqc_report.html
[INFO   ]         multiqc : Data        : WGS_test_multiqc_report_data
[INFO   ]         multiqc : MultiQC complete
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```

This generates a comprehensive [report](https://multiqc.info/) with
a number of metrics summarized across all samples.



Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. multiqc.sh), which uses the input file 'multiqc.in'. For example:



```

#!/bin/bash
module load multiqc || exit 1
multiqc --title 'tenure here i come' /path/to/awsome/project

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=2 --mem=4g multiqc.sh
```







