

document.querySelector('title').textContent = 'fastq\_demux: demultiplexing Illumina FASTQ files based on barcodes in the FASTQ headers';
**fastq\_demux: demultiplexing Illumina FASTQ files based on barcodes in the FASTQ headers.**


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



fastq\_demux is a simple program to demultiplex a FASTQ file or a pair of FASTQ files based on the barcodes present in the FASTQ headers. 




Documentation
* [fastq\_demux GitHub page](https://github.com/Molmed/fastq_demux)


Important Notes
* Module Name: fastq\_demux (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set
	+ **FD\_HOME**  installation directory
	+ **FD\_BIN**       executable directory
	+ **FD\_SRC**       source code directory
	+ **FD\_DATA**  sample data directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive**
[user@cn0911 ~]$**module load fastq\_demux** 
[+] Loading singularity  3.10.5  on cn4208
[+] Loading fastq_demux  20230713
[user@cn0911 ~]$**cp $FD\_DATA\*/\***
[user@cn0911 ~]$**fastq\_demux --R1 dual-index\_Undetermined\_S0\_L001\_I1\_001.fastq.gz --R2 dual-index\_Undetermined\_S0\_L001\_I2\_001.fastq.gz --samplesheet ./samplesheet\_dual\_index.tsv**
known_barcode   count   percent
GGGGGGGG+AGATCTCG       22      44.0
GAAGATTT+TTTACTCT       5       10.0
GAAGATTT+AAAACGCC       3       6.0
unknown_barcode count   percent
GGGGGGGG+AGAACTCG       2       4.0
GGGGGGGG+AGAACGCG       2       4.0
GAAGATTT+CCACTCCG       1       2.0
GGAGATTT+GGGGGGGG       1       2.0
TAAGATTT+TAATCTCT       1       2.0
TAATATTT+TAAACGCT       1       2.0
TTATATTT+TAAACGCT       1       2.0
TCAGGGGG+AGATCTCG       1       2.0
TAATATTT+CCCACGCC       1       2.0
GAAGATTT+AAAACGCG       1       2.0

```
 
End the interactive session:

```

[user@cn0911 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```





