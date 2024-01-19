

document.querySelector('title').textContent = 'tombo on Biowulf';
tombo on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
 |


 Tombo is a suite of tools primarily for the identification of modified nucleotides from nanopore sequencing data. Tombo also provides tools for the analysis and visualization of raw nanopore signal.

 
 Features:
 * Modified Base Detection
	+ Supports both DNA and direct RNA
		- RNA processing details
	+ Three detection algorithms support broad range of applications
		- Alternative model (preferred)
		- Sample comparison
		- De novo
* Reference-anchored raw signal vizualization
* Raw signal analysis python API
* User-friendly model estimation methods with tutorial





### References:


* H. H. Jeong, H. K. Yalamanchili, C. Guo C, J. M. Shulman, and Z. Liu
*De novo Identification of DNA Modifications Enabled by Genome-Guided Nanopore Signal Processing* bioRxiv (2016).
[Journal](http://biorxiv.org/content/early/2017/04/10/094672)


Documentation
* tombo on [GitHub](https://github.com/nanoporetech/tombo)
* tombo on [Documentation](https://nanoporetech.github.io/tombo/)


Important Notes
* Module Name: tombo (see [the modules page](/apps/modules.html)
for more information)
* tombo is multithreaded. Please match the number of threads to the number of allocated CPUs
* Example files in `$TOMBO_TEST_DATA`
* General reference data in `/fdb/igenomes/`



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=5g --cpus-per-task=4 --gres=lscratch:10**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144]$ **cd /lscratch/$SLURM\_JOB\_ID**
[user@cn3144]$ **module load tombo**
[user@cn3144]$ **tar xvfz ${TOMBO\_TEST\_DATA}/${GZFILE} .**
[user@cn3144]$ **export REF="/fdb/igenomes/Saccharomyces\_cerevisiae/UCSC/sacCer3/Sequence/WholeGenomeFasta/genome.fa"**
[user@cn3144]$ **tombo resquiggle $FAST5DIR $REF --processes $SLURM\_CPUS\_PER\_TASK --num-most-common-errors 5 --overwrite**
[18:08:23] Loading minimap2 reference.
[18:08:24] Getting file list.
[18:08:24] Loading default canonical ***** DNA ***** model.
[18:08:24] Re-squiggling reads (raw signal to genomic sequence alignment).
5 most common unsuccessful read types (approx. %):
     2.7% (     93 reads) : Alignment not produced                                                          
     1.0% (     35 reads) : Poor raw to expected signal matching (revert with `tombo filter clear_filters`) 
     0.1% (      5 reads) : Read event to sequence alignment extends beyond bandwidth                       
     -----
     -----
100%|##################################################################| 3621/3621 [05:49<00:00, 10.36it/s]
[18:14:14] Final unsuccessful reads summary (3.7% reads unsuccessfully processed; 133 total reads):
     2.6% (     93 reads) : Alignment not produced                                                          
     1.0% (     35 reads) : Poor raw to expected signal matching (revert with `tombo filter clear_filters`) 
     0.1% (      5 reads) : Read event to sequence alignment extends beyond bandwidth                       
[18:14:14] Saving Tombo reads index to file.

[user@cn3144]$ **tombo detect\_modifications alternative\_model --fast5-basedirs $FAST5DIR \
 --statistics-file-basename native.e\_coli\_sample \
 --alternate-bases dam dcm --processes $SLURM\_CPUS\_PER\_TASK**
    [18:20:03] Parsing Tombo index file(s).
    [18:20:03] Performing alternative model testing.
    [18:20:03] Performing specific alternate base(s) testing.
    [18:20:03] Calculating read coverage regions.
    [18:20:03] Calculating read coverage.
    [18:20:04] Performing modified base detection across genomic regions.
    100%|##############################################################| 2392/2392 [01:10<00:00, 33.89it/s]

[user@cn3144]$ **tombo plot most\_significant --fast5-basedirs $FAST5DIR \
 --statistics-filename native.e\_coli\_sample.dcm.tombo.stats \
 --plot-standard-model --plot-alternate-model dcm \
 --pdf-filename sample.most\_significant\_dcm\_sites.pdf**
    [18:23:05] Loading statistics from file.
    [18:23:06] Parsing Tombo index file(s).
    [18:23:06] Loading default canonical ***** DNA ***** model.
    [18:23:06] Preparing plot data.
    [18:23:08] Plotting.

[user@cn3144]$ **tombo text\_output browser\_files --statistics-filename native.e\_coli\_sample.dam.tombo.stats \
 --file-types dampened\_fraction --browser-file-basename native.e\_coli\_sample.dam**
    [18:23:42] Loading statistics from file.
    [18:23:42] Parsing and outputting statistics wiggles.

[user@cn3144]$ **tombo text\_output browser\_files --fast5-basedirs $FAST5DIR \
 --file-types coverage --browser-file-basename native.e\_coli\_sample**
    [18:24:23] Parsing Tombo index file(s).
    [18:24:23] Getting and writing  coverage bedgraphs.
    [18:24:23] Calculating read coverage regions.
    [18:24:23] Calculating read coverage.   

```



Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch script file (e.g. tombo.sh). For example:



```

#!/bin/bash
cd /lscratch/$SLURM_JOB_ID
module load tombo
tar xvfz ${TOMBO_TEST_DATA}/${GZFILE} .
export REF="/fdb/igenomes/Saccharomyces_cerevisiae/UCSC/sacCer3/Sequence/WholeGenomeFasta/genome.fa"
tombo resquiggle $FAST5DIR $REF --processes $SLURM_CPUS_PER_TASK  --num-most-common-errors 5 --overwrite
....
....

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=4 --mem=10g --gres=lscratch:100 tombo.sh
```







