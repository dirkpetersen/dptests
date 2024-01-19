
















braker:Tool to infer Orthologs from Genome Alignments


  









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







**braker: a pipeline for fully automated prediction of protein coding gene structures with GeneMark-ES/ET and AUGUSTUS in novel eukaryotic genomes**


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



BRAKER3 is the latest pipeline in the BRAKER suite. It enables the usage of RNA-seq and protein data in a fully automated pipeline to train and predict highly reliable genes with GeneMark-ETP and AUGUSTUS. The result of the pipeline is the combined gene set of both gene prediction tools, which only contains genes with very high support from extrinsic evidence.


### References:



 Bruna, T., Hoff, K.J., Lomsadze, A., Stanke, M., & Borodovsky, M. *BRAKER2: Automatic Eukaryotic Genome Annotation with GeneMark-EP+ and AUGUSTUS Supported by a Protein Database* NAR Genomics and Bioinformatics 2021, 3(1):lqaa108, doi: 10.1093/nargab/lqaa108.
   

[PubMed](https://pubmed.ncbi.nlm.nih.gov/33575650/)
[NAR Genomics and Bioinformatics](https://academic.oup.com/nargab/article/3/1/lqaa108/6066535), 2021



Documentation
* [braker GitHub Page](https://github.com/Gaius-Augustus/BRAKER)


Important Notes
* Module Name: braker (see [the modules page](/apps/modules.html) for more information)
* Unusual environment variables set
	+ **BRAKER\_TEST\_DATA** sample data for running braker
	+ BRAKER has to run some steps on a single thread, others can take advantage of multiple threads. But please do not set up more than 8 threads for BRAKER. Otherwise, the rest of threads will be idle.



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive -c 8 --mem=4g --gres=lscratch:10** 
[user@cn3144 ~]$ **module load braker** 
Loading braker  3

[user@cn3144 ]$ **cp -r $BRAKER\_TEST\_DATA/\*.fa .** 

```

run testing data

```

[user@cn3144 ]$ **braker.pl --genome=genome.fa --prot\_seq=proteins.fa --threads $SLURM\_CPUS\_PER\_TASK**


```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. braker.sh). For example:



```

#! /bin/bash

module load braker || exit 1
wget http://bioinf.uni-greifswald.de/augustus/datasets/RNAseq.bam
braker.pl --genome genome.fa --bam RNAseq.bam --threads $SLURM_CPUS_PER_TASK


```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch -c 8 --mem=10g braker.sh
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



  



