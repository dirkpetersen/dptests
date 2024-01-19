

document.querySelector('title').textContent = 'merfin on Biowulf';
merfin on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
 |


From the documentation:

> 
> Improved variant filtering and polishing via k-mer validation
> 





### References:


* G. Formenti et al. *Merfin: improved variant filtering, assembly evaluation and polishing via k-mer validation*.
 Nature Methods 2022. [PubMed](https://pubmed.ncbi.nlm.nih.gov/35361932/) | 
 [PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/pmid/35361932/) | 
 [Journal](https://www.nature.com/articles/s41592-022-01445-y)


Documentation
* merfin on [GitHub](https://github.com/arangrhie/merfin)
* [Best practices](https://github.com/arangrhie/merfin/wiki/Best-practices-for-Merfin)


Important Notes
* Module Name: merfin (see [the modules page](/apps/modules.html) for more information)
* merfin is multithreaded. Match the number of allocated CPUs to the number of threads
* Example files in `$MERFIN_TEST_DATA`



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --cpus-per-task=16 --mem=36g --gres=lscratch:100**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144]$ **cd /lscratch/$SLURM\_JOB\_ID**
[user@cn3144]$ **module load merfin**
[user@cn3144]$ # unpack a meryl k-mer database created from PE 250pb illumina reads from the son of a trio
[user@cn3144]$ # with     meryl count k=21 reads.fastq.gz output HG002.k21.meryl
[user@cn3144]$ # followed by excluding kmers with frequency of 1
[user@cn3144]$ **tar -xzf ${MERFIN\_TEST\_DATA:-none}/HG002.k21.gt1.meryl.tar.gz**
[user@cn3144]$ **cp ${MERFIN\_TEST\_DATA:-none}/{chr20.fasta.gz,ill.vcf.gz} .**
[user@cn3144]$ **merfin -filter -sequence chr20.fasta.gz \
 -memory 34 \
 -threads $SLURM\_CPUS\_PER\_TASK \
 -readmers HG002.k21.gt1.meryl \
 -vcf ill.vcf.gz \
 -output test.merfin**

[user@cn3144]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. merfin.sh), which uses the input file 'merfin.in'. For example:



```

#!/bin/bash
module load merfin
tar -xzf ${MERFIN_TEST_DATA:-none}/HG002.k21.gt1.meryl.tar.gz
cp ${MERFIN_TEST_DATA:-none}/{chr20.fasta.gz,ill.vcf.gz} .
merfin -filter -sequence chr20.fasta.gz  \
    -memory 34 \
    -threads $SLURM_CPUS_PER_TASK \
    -readmers HG002.k21.gt1.meryl \
    -vcf ill.vcf.gz               \
    -output test.merfin

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=16 --mem=36g merfin.sh
```







