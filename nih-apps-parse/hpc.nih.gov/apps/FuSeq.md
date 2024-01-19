

document.querySelector('title').textContent = 'FuSeq: discovering fusion genes from paired-end RNA sequencing data';
**FuSeq: discovering fusion genes from paired-end RNA sequencing data**


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



FuSeq is a fast and accurate method to discover fusion genes based on quasi-mapping to quickly map the reads, extract initial candidates from split reads and fusion equivalence classes of mapped reads, and finally apply multiple filters and statistical tests to get the final candidates. 



### References:



Trung Nghia Vu, Wenjiang Deng, Quang Thinh Trac, Stefano Calza, Woochang Hwang and Yudi Pawitan.   

*A fast detection of fusion genes from paired-end RNA-seq data.*    

[BMC Genomics](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-018-5156-1), 2018, **19**, p.786. <https://doi.org/10.1186/s12864-018-5156-1>.



Documentation
* [FuSeq GitHub Page](https://github.com/nghiavtr/FuSeq)


Important Notes
* Module Name: FuSeq (see [the modules page](/apps/modules.html) for more information)
* Unusual environment variables set
	+ **FUSEQ\_HOME** FuSeq installation directory
	+ **FUSEQ\_BIN** FuSeq executable folder
	+ **FUSEQ\_REF** FuSeq reference directory
	+ **FUSEQ\_DATA** sample data for running FuSeq
	+ **FUSEQ\_R** FuSeq R source files directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=4g --gres=lscratch:10** 
[user@cn3144 ~]$ **module load FuSeq** 
[+] Loading gcc  9.2.0  ...
[+] Loading GSL 2.6 for GCC 9.2.0 ...
[+] Loading openmpi 3.1.4  for GCC 9.2.0
[+] Loading ImageMagick  7.0.8  on cn2367
[+] Loading HDF5  1.10.4
[+] Loading NetCDF 4.7.4_gcc9.2.0
[+] Loading pandoc  2.10.1  on cn2367
[+] Loading pcre2 10.21  ...
[+] Loading R 4.0.0
[+] Loading FuSeq / 1.1.3  ...

```

Create soft links to the sample read data:

```

[user@cn3144 ~]$ **cp -r $FUSEQ\_DATA ./** 

```

Index the sample transcript data if needed, GRCh37\_release75 and GRCh38\_release106 are already indexed at /fdb/fuseq/:

```
 
[user@cn3144 ~]$ **TxIndexer -t Homo\_sapiens.GRCh37.75.cdna.all.fa -o TxIndexer\_idx**
---------------------------------------------------------------

--- FuSeq acknowledges the Sailfish and Rapmap for this indexer ---

---------------------------------------------------------------

writing log to TxIndexer_idx/LogDir/FuSeq_index.log
RapMap Indexer

[Step 1 of 4] : counting k-mers
counted k-mers for 180000 transcriptsElapsed time: 4.19777s

Replaced 0 non-ATCG nucleotides
Clipped poly-A tails from 1401 transcripts
Building rank-select dictionary and saving to disk done
Elapsed time: 0.0257811s
Writing sequence data to file . . . done
Elapsed time: 0.137891s
[info] Building 32-bit suffix array (length of generalized text is 287322183)
Building suffix array . . . success
saving to disk . . . done
Elapsed time: 0.440665s
done
Elapsed time: 41.0055s
processed 287000000 positions
khash had 102501816 keys
saving hash to disk . . . done
Elapsed time: 9.4689s

```

Extract fusion equivalence classes and split reads by running FuSeq with 16 cpus:

```

[user@cn3144 ~]$ **FuSeq -i /fdb/fuseq/TxIndexer\_idx/GRCh37.75/ \
 -l IU -1 SRR064287\_1.fastq -2 SRR064287\_2.fastq -p 16 \
 -g $FUSEQ\_REF/Homo\_sapiens.GRCh37.gtf -o feqDir**

```

- a folder feqDir will be produced.   
   


Finally, discover the fusion genes. This will be done in two steps.   

First, create .sqlite file containing the annotation file:

```

[user@cn3144 ~]$ **Rscript $FUSEQ\_R/createSqlite.R $FUSEQ\_REF/Homo\_sapiens.GRCh37.gtf Homo\_sapiens.GRCh37.sqlite** 

```

Second, use the file Homo\_sapiens.GRCh37.sqlite, together with 
the supporting annotation file Homo\_sapiens.GRCh37.txAnno.RData containing information of paralogs, gene types, etc., to perform the discovery of fusion genes:

```

[user@cn3144 ~]$ **Rscript $FUSEQ\_R/FuSeq.R in=feqDir \
 txfasta=$FUSEQ\_REF/Homo\_sapiens.GRCh37.cdna.all.fa \
 sqlite=Homo\_sapiens.GRCh37.sqlite \
 txanno=Homo\_sapiens.GRCh37.txAnno.RData** 
[user@cn3144 ~]$ **exit**
exit
salloc.exe: Relinquishing job allocation 18320304
[user@biowulf ~]$

```





