
















toga:Tool to infer Orthologs from Genome Alignments


  









|  |  |  |
| --- | --- | --- |
|  | 
[Biowulf High Performance Computing at the NIH](https://hpc.nih.gov)
 | 








[GitHub](https://github.com/NIH-HPC)
[YouTube](https://www.youtube.com/channel/UCx-kNd1kBskYr5KLT9-Erew)
[@nih_hpc](https://twitter.com/nih_hpc)
[RSS Feed](/hpc_RSS.xml)
 |




* [Systems](https://hpc.nih.gov/systems/)
* [Applications](https://hpc.nih.gov/apps/)
* [Reference Data](https://hpc.nih.gov/refdb/)
* [Storage](https://hpc.nih.gov/storage/)
* [User Guides](https://hpc.nih.gov/docs/user_guides.html)
* [Training](https://hpc.nih.gov/training/)
* [User Dashboard](https://hpcnihapps.cit.nih.gov/auth/dashboard/)
* [How To](https://hpc.nih.gov/docs/how_to.html)
* [About](https://hpc.nih.gov/about/)







**toga:Tool to infer Orthologs from Genome Alignments**


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



TOGA is a new method that integrates gene annotation, inferring orthologs and classifying genes as intact or lost.

TOGA implements a novel machine learning based paradigm to infer orthologous genes between related species and to accurately distinguish orthologs from paralogs or processed pseudogenes.

This tutorial explains how to get started using TOGA. It shows how to install and execute TOGA, and how to handle possible issues that may occur.


### References:



Kirilenko BM, Munegowda C, Osipova E, Jebb D, Sharma V, Blumer M, Morales AE, Ahmed AW, Kontopoulos DG, Hilgers L, Lindblad-Toh K, Karlsson EK, Zoonomia Consortium, Hiller M.   

*Integrating gene annotation with orthology inference at scale.*    

[PubMed](https://pubmed.ncbi.nlm.nih.gov/37104600/)
[Science](https://www.science.org/stoken/author-tokens/ST-1161/full), 2023



Documentation
* [toga GitHub Page](https://github.com/hillerlab/TOGA)


Important Notes
* Module Name: toga (see [the modules page](/apps/modules.html) for more information)
* Unusual environment variables set
	+ **TOGA\_GENOME** toga reference directory
	+ **TOGA\_CONFIG** toga config directory
	+ **TOGA\_SUPPLY** toga supplymentary directory
	+ **TOGA\_TEST\_DATA** sample data for running toga



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=4g --gres=lscratch:10** 
[user@cn3144 ~]$ **module load toga** 
[+] Loading java 17.0.3.1  ...
[+] Loading singularity  3.10.5  on cn4292
[+] Loading nextflow 22.10.2
[+] Loading toga  1.1.2

[user@cn3144 ]$ **cp -r $TOGA\_SUPPLY/\* .** 

```

run testing data

```

[user@cn3144 ]$ **toga.py \
 $TOGA\_TEST\_DATA/hg38.mm10.chr11.chain \
 $TOGA\_TEST\_DATA/hg38.genCode27.chr11.bed \
 $TOGA\_GENOME/hg38.2bit \
 $TOGA\_GENOME/mm10.2bit \
 --kt --pn test -i \
 $TOGA\_SUPPLY/hg38.wgEncodeGencodeCompV34.isoforms.txt \
 --nc $TOGA\_CONFIG \
 --cb 3,5 --cjn 500 --u12 supply/hg38.U12sites.tsv --ms \
 --nextflow\_dir ./**


```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. toga.sh). For example:



```

#! /bin/bash

module load toga || exit 1
cp -r $TOGA_SUPPLY/* .

toga.py \
        $TOGA_TEST_DATA/hg38.mm10.chr11.chain \
        $TOGA_TEST_DATA/hg38.genCode27.chr11.bed \
        $TOGA_GENOME/hg38.2bit \
        $TOGA_GENOME/mm10.2bit \
        --kt --pn test -i \
        $TOGA_SUPPLY/hg38.wgEncodeGencodeCompV34.isoforms.txt \
        --nc $TOGA_CONFIG \
        --cb 3,5 --cjn 500 --u12 supply/hg38.U12sites.tsv --ms \
        --nextflow_dir ./


```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch toga.sh
```












[HPC @ NIH](https://hpc.nih.gov)  ~
[Contact](https://hpc.nih.gov/about/contact.html)


[Disclaimer](https://hpc.nih.gov/docs/disclaimer.html) ~ 
[Privacy](https://hpc.nih.gov/docs/privacy.html) ~ 
[Accessibility](https://hpc.nih.gov/docs/accessibility.html) ~ 
[CIT](https://cit.nih.gov/) ~ 
[NIH](https://www.nih.gov/) ~ 
[DHHS](https://www.dhhs.gov/) ~ 
[USA.gov](https://www.firstgov.gov/) ~
[HHS Vulnerability Disclosure](https://www.hhs.gov/vulnerability-disclosure-policy/index.html)



  



