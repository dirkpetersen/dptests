

document.querySelector('title').textContent = "Regtools";
Regtools on Biowulf


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



 RegTools is a set of tools that integrate DNA-seq and RNA-seq data to help interpret mutations in a regulatory and splicing context. 



### References:


* Cotto, K.C., Feng, YY., Ramu, A. et al.
 [**RegTools: Integrated analysis of genomic and transcriptomic data for the discovery of splice-associated variants in cancer**](https://www.nature.com/articles/s41467-023-37266-6)
 *Nat Commun 14, 1589 (2023)*


Documentation
* [Regtools Main Site](https://regtools.readthedocs.io/en/latest/)
* [Regtools Github Site](https://regtools.readthedocs.io/en/latest/)


Important Notes
* Module Name: regtools (see [the modules page](/apps/modules.html) for more information)
* Environment variables set 
	+ REGTOOLS\_HOME* Example files in /usr/local/apps/regtools/1.0.0/regtools/test-data



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive**
[user@cn4338 ~]$ **module load regtools**
[+] Loading regtools 1.0.0  on cn4338 
[+] Loading singularity  3.10.5  on cn4338

```


Example
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Running help command:



```

[user@cn4338] **cp -a /usr/local/apps/regtools/1.0.0/regtools/test-data .**
[user@cn4338 test-data]$ **regtools --help**
Program:        regtools
Version:        1.0.0
Usage:          regtools  [options]
Command:        junctions               Tools that operate on feature junctions (e.g. exon-exon junctions from RNA-seq).
                cis-ase                 Tools related to allele specific expression in cis.
                cis-splice-effects      Tools related to splicing effects of variants.
                variants                Tools that operate on variants.


```

Create a batch input file (e.g. regtools.sh). For example:



```

#!/bin/bash
module load regtools

regtools cis-splice-effects identify -s FR variants.vcf alignments.bam ref.fa annotations.gtf

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] regtools.sh
```

cis-splice-effects:



```

[user@cn4338 test-data]$**regtools cis-splice-effects identify -s RF \
 test4.vcf.gz \
 test\_hcc1395.bam \
 test\_chr22.fa \
 test\_ensemble\_chr22.gtf**

  Program:        regtools
  Version:        1.0.0
  Variant file: test4.vcf.gz
  Alignment file: test_hcc1395.bam
  Reference fasta file: test_chr22.fa
  Annotation file: test_ensemble_chr22.gtf
  
  exonic_min_distance_ is 3
  
  chrom   start   end     name    score   strand  splice_site     acceptors_skipped       exons_skipped   donors_skipped  anchor  known_donor     known_acceptor  known_junction  gene_names      gene_ids transcripts      variant_info


```



|  |  |
| --- | --- |
|  For more examples, please visit the [Regtools Documentation Page](https://regtools.readthedocs.io/en/latest/) | |








