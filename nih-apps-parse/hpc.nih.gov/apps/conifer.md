

document.querySelector('title').textContent = 'conifer: finding copy number variants and genotyping the copy-number of duplicated genes.';
conifer: finding copy number variants and genotyping the copy-number of duplicated genes.


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



CoNIFER (copy number inference from exome reads) uses exome sequencing data to find copy number variants (CNVs) and genotype the copy-number of duplicated genes. It can reliably be
used to discover disruptive genic CNVs missed by standard approaches and should have broad
application in human genetic studies of disease.
  

  




### References:


* Niklas Krumm, Peter H. Sudmant, Arthur Ko, Brian J. O’Roak, Maika Malig, Bradley P. Coe, NHLBI Exome Sequencing Project, Aaron R. Quinlan, Deborah A. Nickerson, Evan E. Eichler   

*Copy number variation detection and genotyping from exome sequence data*    

[Genome Res. 2012. 22: 1525-1532](https://genome.cshlp.org/content/early/2012/05/14/gr.138115.112.full.pdf+html)


Documentation
* [CoNIFER Sourceforge page](https://sourceforge.net/projects/conifer/)
* [CoNIFER tutorial](https://conifer.sourceforge.net/tutorial.html)


Important Notes
* Module Name: conifer (see [the modules page](/apps/modules.html) for more information)
  
* Unusual environment variables set
	+ **CONIFER\_HOME**  installation directory
	+ **CONIFER\_SRC**  source code directory
	+ **CONIFER\_DATA**  sample data directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive**
[user@cn3144 ~]$ **module load conider**
[+] Loading singularity  3.10.5  on cn4198
[+] Loading conifer  0.2.2

```

  

The basdic usage of the conifer application is as follows:

```

[user@cn3144 ~]$ **python conifer.py -h** 
usage: CoNIFER [-h] [--version] {rpkm,analyze,export,plot,call,plotcalls} ...

This is CoNIFER 0.2.2 (Copy Number Inference From Exome Reads), designed to
detect and genotype CNVs and CNPs from exome sequence read-depth data. See
Krumm et al., Genome Research (2012) doi:10.1101/gr.138115.112 Niklas Krumm,
2012 nkrumm@uw.edu

positional arguments:
  {rpkm,analyze,export,plot,call,plotcalls}
                        Command to be run.
    rpkm                Create an RPKM file from a BAM file
    analyze             Basic CoNIFER analysis. Reads a directory of RPKM
                        files and a probe list and outputs a HDF5 file
                        containing SVD-ZRPKM values.
    export              Export SVD-ZRPKM values from a HDF5 file to bed or vcf
                        format.
    plot                Plot SVD-ZRPKM values using matplotlib
    call                Very rudimentary caller for CNVs using SVD-ZRPKM
                        thresholding.
    plotcalls           Make basic plots from call file from "call" command.

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
[user@cn3144 ~]$ **cp $CONIFER\_DATA/\* .**
[user@cn3144 ~]$ **conifer.py analyze --probes probes.txt --rpkm\_dir RPKM\_data/ --output analysis.hdf5 --svd 6 --write\_svals singular\_values.txt**
[INIT] Successfully read in 194080 probes from probes.txt
[INIT] Found 26 RPKM files in RPKM_data/
[INIT] Mapping file to sampleID: RPKM_data/NA12878.txt --> NA12878
[INIT] Mapping file to sampleID: RPKM_data/NA15510.txt --> NA15510
[INIT] Mapping file to sampleID: RPKM_data/NA18507.txt --> NA18507
[INIT] Mapping file to sampleID: RPKM_data/NA18517.txt --> NA18517
[INIT] Mapping file to sampleID: RPKM_data/NA18555.txt --> NA18555
[INIT] Mapping file to sampleID: RPKM_data/NA18956.txt --> NA18956
[INIT] Mapping file to sampleID: RPKM_data/NA19129.txt --> NA19129
[INIT] Mapping file to sampleID: RPKM_data/NA19240.txt --> NA19240
[INIT] Mapping file to sampleID: RPKM_data/Trio1.fa.txt --> Trio1.fa
[INIT] Mapping file to sampleID: RPKM_data/Trio1.mo.txt --> Trio1.mo
[INIT] Mapping file to sampleID: RPKM_data/Trio1.p1.txt --> Trio1.p1
[INIT] Mapping file to sampleID: RPKM_data/Trio2.fa.txt --> Trio2.fa
[INIT] Mapping file to sampleID: RPKM_data/Trio2.mo.txt --> Trio2.mo
[INIT] Mapping file to sampleID: RPKM_data/Trio2.p1.txt --> Trio2.p1
[INIT] Mapping file to sampleID: RPKM_data/Trio3.fa.txt --> Trio3.fa
[INIT] Mapping file to sampleID: RPKM_data/Trio3.mo.txt --> Trio3.mo
[INIT] Mapping file to sampleID: RPKM_data/Trio3.p1.txt --> Trio3.p1
[INIT] Mapping file to sampleID: RPKM_data/Trio4.fa.txt --> Trio4.fa
[INIT] Mapping file to sampleID: RPKM_data/Trio4.mo.txt --> Trio4.mo
[INIT] Mapping file to sampleID: RPKM_data/Trio4.p1.txt --> Trio4.p1
[INIT] Mapping file to sampleID: RPKM_data/Trio5.fa.txt --> Trio5.fa
[INIT] Mapping file to sampleID: RPKM_data/Trio5.mo.txt --> Trio5.mo
[INIT] Mapping file to sampleID: RPKM_data/Trio5.p1.txt --> Trio5.p1
[INIT] Mapping file to sampleID: RPKM_data/Trio6.fa.txt --> Trio6.fa
[INIT] Mapping file to sampleID: RPKM_data/Trio6.mo.txt --> Trio6.mo
[INIT] Mapping file to sampleID: RPKM_data/Trio6.p1.txt --> Trio6.p1
[INIT] Successfully read RPKM data for sampleID: Trio6.mo
[INIT] Successfully read RPKM data for sampleID: Trio4.mo
[INIT] Successfully read RPKM data for sampleID: Trio3.fa
[INIT] Successfully read RPKM data for sampleID: Trio1.fa
[INIT] Successfully read RPKM data for sampleID: Trio2.mo
[INIT] Successfully read RPKM data for sampleID: Trio3.p1
[INIT] Successfully read RPKM data for sampleID: Trio4.fa
[INIT] Successfully read RPKM data for sampleID: Trio5.fa
[INIT] Successfully read RPKM data for sampleID: Trio2.fa
[INIT] Successfully read RPKM data for sampleID: Trio1.mo
[INIT] Successfully read RPKM data for sampleID: NA18517
[INIT] Successfully read RPKM data for sampleID: NA12878
[INIT] Successfully read RPKM data for sampleID: Trio5.p1
[INIT] Successfully read RPKM data for sampleID: NA18507
[INIT] Successfully read RPKM data for sampleID: NA18555
[INIT] Successfully read RPKM data for sampleID: NA15510
[INIT] Successfully read RPKM data for sampleID: NA18956
[INIT] Successfully read RPKM data for sampleID: Trio2.p1
[INIT] Successfully read RPKM data for sampleID: Trio3.mo
[INIT] Successfully read RPKM data for sampleID: NA19240
[INIT] Successfully read RPKM data for sampleID: Trio6.fa
[INIT] Successfully read RPKM data for sampleID: Trio5.mo
[INIT] Successfully read RPKM data for sampleID: Trio4.p1
[INIT] Successfully read RPKM data for sampleID: Trio1.p1
[INIT] Successfully read RPKM data for sampleID: NA19129
[INIT] Successfully read RPKM data for sampleID: Trio6.p1
[INIT] Finished reading RPKM files. Total number of samples in experiment: 26 (0 failed to read properly)
[INIT] Attempting to process chromosomes:  chr1, chr2, chr3, chr4, chr5, chr6, chr7, chr8, chr9, chr10, chr11, chr12, chr13, chr14, chr15, chr16, chr17, chr18, chr19, chr20, chr21, chr22, chrX, chrY
[RUNNING: chr1] Now on: chr1
[RUNNING: chr1] Found 19822 probes; probeID range is [0-19822]
[RUNNING: chr1] Calculating median RPKM
[RUNNING: chr1] Masking 412 probes with median RPKM < 1.000000
[RUNNING: chr1] Calculating ZRPKM scores...
[RUNNING: chr1] SVD decomposition...
[RUNNING: chr1] Saving SVD-ZRPKM values
...
[RUNNING: chr2] Now on: chr2
[RUNNING: chr2] Found 14782 probes; probeID range is [19822-34604]
[RUNNING: chr2] Calculating median RPKM
[RUNNING: chr2] Masking 280 probes with median RPKM < 1.000000
[RUNNING: chr2] Calculating ZRPKM scores...
[RUNNING: chr2] SVD decomposition...
[RUNNING: chr2] Saving SVD-ZRPKM values
[RUNNING: chr3] Now on: chr3
[RUNNING: chr3] Found 11714 probes; probeID range is [34604-46318]
[RUNNING: chr3] Calculating median RPKM
[RUNNING: chr3] Masking 207 probes with median RPKM < 1.000000
[RUNNING: chr3] Calculating ZRPKM scores...
[RUNNING: chr3] SVD decomposition...
[RUNNING: chr3] Saving SVD-ZRPKM values
[RUNNING: chr4] Now on: chr4
[RUNNING: chr4] Found 7549 probes; probeID range is [46318-53867]
[RUNNING: chr4] Calculating median RPKM
[RUNNING: chr4] Masking 142 probes with median RPKM < 1.000000
[RUNNING: chr4] Calculating ZRPKM scores...
[RUNNING: chr4] SVD decomposition...
[RUNNING: chr4] Saving SVD-ZRPKM values
[RUNNING: chr5] Now on: chr5
[RUNNING: chr5] Found 8602 probes; probeID range is [53867-62469]
[RUNNING: chr5] Calculating median RPKM
[RUNNING: chr5] Masking 168 probes with median RPKM < 1.000000
[RUNNING: chr5] Calculating ZRPKM scores...
[RUNNING: chr5] SVD decomposition...
[RUNNING: chr5] Saving SVD-ZRPKM values
[RUNNING: chr6] Now on: chr6
[RUNNING: chr6] Found 9221 probes; probeID range is [62469-71690]
[RUNNING: chr6] Calculating median RPKM
[RUNNING: chr6] Masking 284 probes with median RPKM < 1.000000
[RUNNING: chr6] Calculating ZRPKM scores...
[RUNNING: chr6] SVD decomposition...
[RUNNING: chr6] Saving SVD-ZRPKM values
[RUNNING: chr7] Now on: chr7
[RUNNING: chr7] Found 9213 probes; probeID range is [71690-80903]
[RUNNING: chr7] Calculating median RPKM
[RUNNING: chr7] Masking 217 probes with median RPKM < 1.000000
[RUNNING: chr7] Calculating ZRPKM scores...
[RUNNING: chr7] SVD decomposition...
[RUNNING: chr7] Saving SVD-ZRPKM values
[RUNNING: chr8] Now on: chr8
[RUNNING: chr8] Found 6597 probes; probeID range is [80903-87500]
[RUNNING: chr8] Calculating median RPKM
[RUNNING: chr8] Masking 137 probes with median RPKM < 1.000000
[RUNNING: chr8] Calculating ZRPKM scores...
[RUNNING: chr8] SVD decomposition...
[RUNNING: chr8] Saving SVD-ZRPKM values
[RUNNING: chr9] Now on: chr9
[RUNNING: chr9] Found 7887 probes; probeID range is [87500-95387]
[RUNNING: chr9] Calculating median RPKM
[RUNNING: chr9] Masking 216 probes with median RPKM < 1.000000
[RUNNING: chr9] Calculating ZRPKM scores...
[RUNNING: chr9] SVD decomposition...
[RUNNING: chr9] Saving SVD-ZRPKM values
[RUNNING: chr10] Now on: chr10
[RUNNING: chr10] Found 8018 probes; probeID range is [95387-103405]
[RUNNING: chr10] Calculating median RPKM
[RUNNING: chr10] Masking 141 probes with median RPKM < 1.000000
[RUNNING: chr10] Calculating ZRPKM scores...
[RUNNING: chr10] SVD decomposition...
[RUNNING: chr10] Saving SVD-ZRPKM values
[RUNNING: chr11] Now on: chr11
[RUNNING: chr11] Found 10857 probes; probeID range is [103405-114262]
[RUNNING: chr11] Calculating median RPKM
[RUNNING: chr11] Masking 194 probes with median RPKM < 1.000000
[RUNNING: chr11] Calculating ZRPKM scores...
[RUNNING: chr11] SVD decomposition...
[RUNNING: chr11] Saving SVD-ZRPKM values
[RUNNING: chr12] Now on: chr12
[RUNNING: chr12] Found 10976 probes; probeID range is [114262-125238]
[RUNNING: chr12] Calculating median RPKM
[RUNNING: chr12] Masking 172 probes with median RPKM < 1.000000
[RUNNING: chr12] Calculating ZRPKM scores...
[RUNNING: chr12] SVD decomposition...
[RUNNING: chr12] Saving SVD-ZRPKM values
[RUNNING: chr13] Now on: chr13
[RUNNING: chr13] Found 3475 probes; probeID range is [125238-128713]
[RUNNING: chr13] Calculating median RPKM
[RUNNING: chr13] Masking 81 probes with median RPKM < 1.000000
[RUNNING: chr13] Calculating ZRPKM scores...
[RUNNING: chr13] SVD decomposition...
[RUNNING: chr13] Saving SVD-ZRPKM values
[RUNNING: chr14] Now on: chr14
[RUNNING: chr14] Found 6098 probes; probeID range is [128713-134811]
[RUNNING: chr14] Calculating median RPKM
[RUNNING: chr14] Masking 113 probes with median RPKM < 1.000000
[RUNNING: chr14] Calculating ZRPKM scores...
[RUNNING: chr14] SVD decomposition...
[RUNNING: chr14] Saving SVD-ZRPKM values
[RUNNING: chr15] Now on: chr15
[RUNNING: chr15] Found 6766 probes; probeID range is [134811-141577]
[RUNNING: chr15] Calculating median RPKM
[RUNNING: chr15] Masking 124 probes with median RPKM < 1.000000
[RUNNING: chr15] Calculating ZRPKM scores...
[RUNNING: chr15] SVD decomposition...
[RUNNING: chr15] Saving SVD-ZRPKM values
[RUNNING: chr16] Now on: chr16
[RUNNING: chr16] Found 8262 probes; probeID range is [141577-149839]
[RUNNING: chr16] Calculating median RPKM
[RUNNING: chr16] Masking 225 probes with median RPKM < 1.000000
[RUNNING: chr16] Calculating ZRPKM scores...
[RUNNING: chr16] SVD decomposition...
[RUNNING: chr16] Saving SVD-ZRPKM values
[RUNNING: chr17] Now on: chr17
[RUNNING: chr17] Found 11596 probes; probeID range is [149839-161435]
[RUNNING: chr17] Calculating median RPKM
[RUNNING: chr17] Masking 238 probes with median RPKM < 1.000000
[RUNNING: chr17] Calculating ZRPKM scores...
[RUNNING: chr17] SVD decomposition...
[RUNNING: chr17] Saving SVD-ZRPKM values
[RUNNING: chr18] Now on: chr18
[RUNNING: chr18] Found 2988 probes; probeID range is [161435-164423]
[RUNNING: chr18] Calculating median RPKM
[RUNNING: chr18] Masking 47 probes with median RPKM < 1.000000
[RUNNING: chr18] Calculating ZRPKM scores...
[RUNNING: chr18] SVD decomposition...
[RUNNING: chr18] Saving SVD-ZRPKM values
[RUNNING: chr19] Now on: chr19
[RUNNING: chr19] Found 11162 probes; probeID range is [164423-175585]
[RUNNING: chr19] Calculating median RPKM
[RUNNING: chr19] Masking 338 probes with median RPKM < 1.000000
[RUNNING: chr19] Calculating ZRPKM scores...
[RUNNING: chr19] SVD decomposition...
[RUNNING: chr19] Saving SVD-ZRPKM values
[RUNNING: chr20] Now on: chr20
[RUNNING: chr20] Found 4931 probes; probeID range is [175585-180516]
[RUNNING: chr20] Calculating median RPKM
[RUNNING: chr20] Masking 166 probes with median RPKM < 1.000000
[RUNNING: chr20] Calculating ZRPKM scores...
[RUNNING: chr20] SVD decomposition...
[RUNNING: chr20] Saving SVD-ZRPKM values
[RUNNING: chr21] Now on: chr21
[RUNNING: chr21] Found 2066 probes; probeID range is [180516-182582]
[RUNNING: chr21] Calculating median RPKM
[RUNNING: chr21] Masking 71 probes with median RPKM < 1.000000
[RUNNING: chr21] Calculating ZRPKM scores...
[RUNNING: chr21] SVD decomposition...
[RUNNING: chr21] Saving SVD-ZRPKM values
[RUNNING: chr22] Now on: chr22
[RUNNING: chr22] Found 4159 probes; probeID range is [182582-186741]
[RUNNING: chr22] Calculating median RPKM
[RUNNING: chr22] Masking 134 probes with median RPKM < 1.000000
[RUNNING: chr22] Calculating ZRPKM scores...
[RUNNING: chr22] SVD decomposition...
[RUNNING: chr22] Saving SVD-ZRPKM values
[RUNNING: chr23] Now on: chrX
[RUNNING: chr23] Found 6874 probes; probeID range is [186741-193615]
[RUNNING: chr23] Calculating median RPKM
[RUNNING: chr23] Masking 156 probes with median RPKM < 1.000000
[RUNNING: chr23] Calculating ZRPKM scores...
[RUNNING: chr23] SVD decomposition...
[RUNNING: chr23] Saving SVD-ZRPKM values
[RUNNING: chr24] Now on: chrY
[RUNNING: chr24] Found 465 probes; probeID range is [193615-194080]
[RUNNING: chr24] Calculating median RPKM
[RUNNING: chr24] Masking 158 probes with median RPKM < 1.000000
[RUNNING: chr24] Calculating ZRPKM scores...
[RUNNING: chr24] SVD decomposition...
[RUNNING: chr24] Saving SVD-ZRPKM values
[RUNNING] Saving sampleIDs to file...
[FINISHED]

```

  

End the interactive session:

```

[user@cn3111 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```





