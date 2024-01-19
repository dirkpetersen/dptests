

document.querySelector('title').textContent = 'SRST2: Short Read Sequence Typing for Bacterial Pathogens ';
**SRST2: Short Read Sequence Typing for Bacterial Pathogens** 


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



SRST2 is a a read mapping-based tool for rapid molecular typing of bacterial pathogens. 
It allows fast and accurate detection of genes, alleles and multi-locus
sequence types (MLST) from WGS data. SRST2 is highly accurate and outperforms assembly-based methods in terms of both gene detection and allele assignment.




### References:


* M.Inouye, H.Dashnow, L.-A.Raven, M.B.Schultz, B.J.Pope, 
 T.Tomita, J.Zobel and K.E.Holt   

*SRST2: Rapid genomic surveillance for public health and hospital microbiology labs*  

[Genome Medicine](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-014-0090-6)  2014, **6:** 90


Documentation
* [SRST2 GitHub page](https://github.com/katholt/srst2)


Important Notes
* Module Name: SRST2 (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set
	+ **SRST2\_HOME**  installation directory
	+ **SRST2\_BIN**       executable directory
	+ **SRST2\_SRC**       source code directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive** 
[user@cn3200 ~]$ **module load srst2** 
[+] Loading bowtie  2-2.2.6
[+] Loading srst2 0.2.0  ...

```

Run a sample command:

```

[user@cn3200 ~]$ **getmlst.py --species "Staphylococcus aureus"**

  For SRST2, remember to check what separator is being used in this allele database

  Looks like --mlst_delimiter '_'

  >arcC_1  --> -->   ('arcC', '_', '1')

  Suggested srst2 command for use with this MLST database:

    srst2 --output test --input_pe *.fastq.gz --mlst_db Staphylococcus_aureus.fasta --mlst_definitions saureus.txt --mlst_delimiter '_'

```

The following files will be produced in the current folder:

```

[user@cn3200 ~]$ **tree .**
.
|-- Staphylococcus_aureus.fasta
|-- arcC.tfa
|-- aroE.tfa
|-- glpF.tfa
|-- gmk.tfa
|-- mlst_data_download_Staphylococcus_aureus_2019-08-20.log
|-- pta.tfa
|-- saureus.txt
|-- tpi.tfa
`-- yqiL.tfa

0 directories, 10 files

```

Run another sample command:

```

[user@cn3200 ~]$ **getmlst.py --species "Staphylococcus epidermidis" --repository\_url http://pubmlst.org/data/dbases.xml**

  For SRST2, remember to check what separator is being used in this allele database

  Looks like --mlst_delimiter '_'

  >arcC_1  --> -->   ('arcC', '_', '1')

  Suggested srst2 command for use with this MLST database:

    srst2 --output test --input_pe *.fastq.gz --mlst_db Staphylococcus_epidermidis.fasta --mlst_definitions sepidermidis.txt --mlst_delimiter '_'

```

The output files are as follows:

```

[user@cn3200 ~]$ **tree .**
.
|-- Staphylococcus_epidermidis.fasta
|-- arcC.tfa
|-- aroE.tfa
|-- gtr.tfa
|-- mlst_data_download_Staphylococcus_epidermidis_2019-08-20.log
|-- mutS.tfa
|-- pyrR.tfa
|-- sepidermidis.txt
|-- tpiA.tfa
`-- yqiL.tfa

0 directories, 10 files

```

End the interactive session:

```

[user@cn3200 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```





