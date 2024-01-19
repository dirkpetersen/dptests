

document.querySelector('title').textContent = 'mosaicforecast on Biowulf';
mosaicforecast on Biowulf


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



MosaicForecast is a machine learning method that leverages read-based phasing and read-level features to accurately detect mosaic SNVs (SNPs, small indels) from NGS data. It builds on existing algorithms to achieve a multifold increase in specificity. 








### References:


* Dou Y, Kwon M, Rodin RE, Cortés-Ciriano I, Doan R, Luquette LJ, Galor A, Bohrson C, Walsh CA, Park PJ. *Accurate detection of mosaic variants in sequencing data without matched controls* Nat Biotechnol. 2020 Mar;38(3):314-319. doi: 10.1038/s41587-019-0368-8. Epub 2020 Jan 6.
 [PubMed](https://www.ncbi.nlm.nih.gov/pubmed/31907404) | 
 [Journal](https://www.nature.com/articles/s41587-019-0368-8)


Documentation
* mosaicforecast Main Site:[Main Site](https://github.com/parklab/MosaicForecast)


Important Notes
* Module Name: mosaicforecast (see [the modules page](/apps/modules.html) for more information)
 * Mosaicforecast was installed as a container. There are 6 functions included in the packages, and please call them directly (without python or Rscript).
 
```

        Phase.py
	ReadLevel_Features_extraction.py
	Prediction.R
	Train_RFmodel.R
	PhasingRefine.R
	MuTect2-PoN_filter.py
	
```
* Environment variables set 
	+ $MOSAIC\_TESTDATA #include demo and umap\_mappability(bigWig file,k=24) from hg19
	 + $MOSAIC\_MODEL #include pre-trained model to predic genotypes, please choose the one match your sequence depths.



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program.   
Sample session (user input in **bold**):



```

[user@biowulf]$ **sinteractive --cpus-per-task=2 --mem=2G**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load mosaicforecast**
[user@cn3144 ~]$ **Phase.py**
Usage: python Phase.py bam_dir output_dir ref_fasta input_positions(file format:chr pos-1 pos ref alt sample, sep=\t) min_dp_inforSNPs(int) Umap_mappability(bigWig file,k=24) n_threads_parallel sequencing_file_format(bam/cram)

Note:
1. Name of bam files should be "sample.bam" under the bam_dir, and there should be corresponding index files.
2. There should be a fai file under the same dir of the fasta file (samtools faidx input.fa).
3. The "min_dp_inforSNPs" is the minimum depth of coverage of trustworthy neaby het SNPs.
4. Bam file is preferred than cram file, as the program would run much more slowly if using cram format.

[user@cn3144 ~]$ **mkdir mosaicforecast\_test && cd mosaicforecast**
[user@cn3144 ~]$ **cp -r ${MOSAIC\_TESTDATA:-none}/\* .**
[user@cn3144 ~]$  **Phase.py ./demo/ test\_out \
 /fdb/GATK\_resource\_bundle/b37-2.8/human\_g1k\_v37\_decoy.fasta \
 ./demo/test.input 20 k24.umap.wg.bw 2 bam**
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. mosaicforecast.sh). For example:



```

#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G
#SBATCH --time=2:00:00
#SBATCH --partition=norm

set -e
module load mosaicforecast
cp -r ${MOSAIC_TESTDATA:-none}/* .
cp -r ${MOSAIC_MODEL:-none}/* .
Prediction.R demo/test.SNP.features models_trained/250xRFmodel_addRMSK_Refine.rds Refine test.SNP.predictions

```

 Submit the job:

```
sbatch mosaicforecast.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. job.swarm). For example:



```


       Prediction.R demo/test.SNP.features models_trained/250xRFmodel_addRMSK_Refine.rds Refine SNP.predictions   
       Prediction.R demo/test.DEL.features models_trained/deletions_250x.RF.rds Phase DEL.predictions

    
```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f job.swarm [-g #] --module mosaicforecast
```

where
 

|  |  |
| --- | --- |
| -g *#*  | Number of Gigabytes of memory required for each process
 (1 line in the swarm command file)  |
| --module  | Loads the module for each subjob in the swarm  |










