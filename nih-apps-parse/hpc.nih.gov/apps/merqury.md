

document.querySelector('title').textContent = 'merqury on Biowulf';
merqury on Biowulf


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



Description from the GitHub repo




> 
>  Evaluate genome assemblies with k-mers and more
> 
>  Often, genome assembly projects have illumina whole genome sequencing reads
>  available for the assembled individual.
> 
>  The k-mer spectrum of this read set can be used for independently evaluating
>  assembly quality without the need of a high quality reference.
> 
>  Merqury provides a set of tools for this purpose.
> 


### References:


* A. Rhie, B. P. Walenz, S. Koren et al. *Merqury: reference-free quality, completeness, and phasing assessment for genome assemblies*. Genome Biol (2020). [PubMed](https://pubmed.ncbi.nlm.nih.gov/32928274/) | 
 [PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/pmid/32928274/) | 
 [Journal](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02134-9)


Documentation
* merqury on [GitHub](https://github.com/marbl/merqury)


Important Notes
* Module Name: merqury (see [the modules page](/apps/modules.html) for more information)
* Best run with approximately 6 CPUs in local mode.
* Example files in `$MERQURY_TEST_DATA`



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=10g --cpus-per-task=6 --gres=lscratch:20**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144]$ **cd /lscratch/$SLURM\_JOB\_ID**
[user@cn3144]$ **module load merqury**
[user@cn3144]$ **cp ${MERQURY\_TEST\_DATA:-none}/\* .**
[user@cn3144]$ **for f in \*.tar.gz ; do tar -xzf $f ; done**
[user@cn3144]$ **merqury.sh F1.k18.meryl col0.hapmer.meryl cvi0.hapmer.meryl athal\_COL.fasta athal\_CVI.fasta test**
[user@cn3144]$ **ls test\***
test.athal_COL.100_20000.phased_block.bed       test.athal_COL.spectra-cn.fl.png                test.athal_CVI.cvi0.hapmer.wig
test.athal_COL.100_20000.phased_block.blob.png  test.athal_COL.spectra-cn.hist                  test.athal_CVI.only.hist
test.athal_COL.100_20000.phased_block.counts    test.athal_COL.spectra-cn.ln.png                test.athal_CVI.qv
test.athal_COL.100_20000.phased_block.sizes     test.athal_COL.spectra-cn.st.png                test.athal_CVI.sort.bed
test.athal_COL.100_20000.phased_block.stats     test.athal_CVI.100_20000.phased_block.bed       test.athal_CVI.spectra-cn.fl.png
test.athal_COL.100_20000.switch.bed             test.athal_CVI.100_20000.phased_block.blob.png  test.athal_CVI.spectra-cn.hist
test.athal_COL.100_20000.switches.txt           test.athal_CVI.100_20000.phased_block.counts    test.athal_CVI.spectra-cn.ln.png
test.athal_COL.block.N.png                      test.athal_CVI.100_20000.phased_block.sizes     test.athal_CVI.spectra-cn.st.png
test.athal_COL.col0.hapmer.spectra-cn.fl.png    test.athal_CVI.100_20000.phased_block.stats     test.completeness.stats
test.athal_COL.col0.hapmer.spectra-cn.ln.png    test.athal_CVI.100_20000.switch.bed             test.dist_only.hist
test.athal_COL.col0.hapmer.spectra-cn.st.png    test.athal_CVI.100_20000.switches.txt           test.hapmers.blob.png
test.athal_COL.col0.hapmer.spectra-hap-cn.hist  test.athal_CVI.block.N.png                      test.hapmers.count
test.athal_COL.col0.hapmer.wig                  test.athal_CVI.col0.hapmer.spectra-cn.fl.png    test.only.hist
test.athal_COL.contig.sizes                     test.athal_CVI.col0.hapmer.spectra-cn.ln.png    test.qv
test.athal_COL.continuity.N.png                 test.athal_CVI.col0.hapmer.spectra-cn.st.png    test.spectra-asm.fl.png
test.athal_COL.cvi0.hapmer.spectra-cn.fl.png    test.athal_CVI.col0.hapmer.spectra-hap-cn.hist  test.spectra-asm.hist
test.athal_COL.cvi0.hapmer.spectra-cn.ln.png    test.athal_CVI.col0.hapmer.wig                  test.spectra-asm.ln.png
test.athal_COL.cvi0.hapmer.spectra-cn.st.png    test.athal_CVI.contig.sizes                     test.spectra-asm.st.png
test.athal_COL.cvi0.hapmer.spectra-hap-cn.hist  test.athal_CVI.continuity.N.png                 test.spectra-cn.fl.png
test.athal_COL.cvi0.hapmer.wig                  test.athal_CVI.cvi0.hapmer.spectra-cn.fl.png    test.spectra-cn.hist
test.athal_COL.only.hist                        test.athal_CVI.cvi0.hapmer.spectra-cn.ln.png    test.spectra-cn.ln.png
test.athal_COL.qv                               test.athal_CVI.cvi0.hapmer.spectra-cn.st.png    test.spectra-cn.st.png
test.athal_COL.sort.bed                         test.athal_CVI.cvi0.hapmer.spectra-hap-cn.hist
[user@cn3144]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. merqury.sh), which uses the test data in `$MERQURY_TEST_DATA`. For example:



```

#!/bin/bash
module load merqury/1.3
cp ${MERQURY_TEST_DATA:-none}/* .
for f in *.tar.gz ; do tar -xzf $f ; done
merqury.sh F1.k18.meryl col0.hapmer.meryl cvi0.hapmer.meryl athal_COL.fasta athal_CVI.fasta test

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=6 --mem=10g merqury.sh
```







