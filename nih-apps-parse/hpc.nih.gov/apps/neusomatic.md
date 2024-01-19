

document.querySelector('title').textContent = 'NeuSomatic on Biowulf';
NeuSomatic on Biowulf


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



NeuSomatic is based on deep convolutional neural networks for accurate somatic mutation detection. With properly trained models, it can robustly perform across sequencing platforms, strategies, and conditions. NeuSomatic summarizes and augments sequence alignments in a novel way and incorporates multi-dimensional features to capture variant signals effectively. It is not only a universal but also accurate somatic mutation detection method.



### References:


* [Sahraeian, Sayed Mohammad Ebrahim, et al. "Deep convolutional neural networks for accurate somatic mutation detection." *Nature communications* 10.1 (2019): 1041.](https://www.nature.com/articles/s41467-019-09027-x)
* [Sahraeian, Sayed Mohammad Ebrahim, et al. "Robust Cancer Mutation Detection with Deep Learning Models Derived from Tumor-Normal Sequencing Data." *BioRxiv* (2019): 667261.](https://www.biorxiv.org/content/10.1101/667261v2)


Documentation
* [Neusomatic Main Site](https://github.com/bioinform/neusomatic)


Important Notes
* Module Name: neusomatic (see [the modules page](/apps/modules.html) for more information)
* This is a GPU app. (Will default to CPU if $CUDA\_VISIBLE\_DEVICES is unset or empty.)
* Newer versions of neusomatic are installed as a container. 
* Environment variables set when the module is loaded. (NOTE: these refer to locations within the neusomatic container and do not exist on the host. They can be passed as arguments to CLI options at runtime. See examples below.)
	+ NEUSOMATIC\_BIN=/opt/neusomatic/neusomatic/bin
	+ NEUSOMATIC\_MODELS=/opt/neusomatic/neusomatic/models
	+ NEUSOMATIC\_SCAN\_ALIGNMENTS=/opt/neusomatic/neusomatic/bin/scan\_alignments* You can use the --checkpoint flag to start training neusomatic from a pre-trained model. See the main site above for more information about this and other options. 
* Example files are maintained in the [GitHub repository under "tests".](https://github.com/bioinform/neusomatic/tree/master/test) A modified version of the run\_test.sh script that is better attuned to the Biowulf environment can be found at /usr/local/apps/neusomatic/0.2.1/neusomatic/test* You must ensure that the name of the reference matches that in the header of the \*.bam file or you will receive difficult to understand errors.
* The neusomatic directory can be copied outside of your home directory. 
* You can run your own data by replacing the input files, adjusting the parameters and choosing the proper pretrained model in the example script. This may be the easiest way to get started analyzing your own data. 
* For best performance, you probably need to train a model on your own data set. 
* If you receive errors like assert(fasta\_file.fetch((c), p - 1, p - 1 + len(r)).upper() == r) you should check to make sure that your truth\_vcf file does not container lowercase letters in the REF and ALT fields.



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --gres=gpu:k80:1,lscratch:10 --mem=20g -c14**
salloc.exe: Pending job allocation 40260664
salloc.exe: job 40260664 queued and waiting for resources
salloc.exe: job 40260664 has been allocated resources
salloc.exe: Granted job allocation 40260664
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn4192 are ready for job
srun: error: x11: no local DISPLAY defined, skipping

[user@cn4192 ~]$ **cp -r /usr/local/apps/neusomatic/0.2.1/neusomatic ~/**

[user@cn4192 ~]$ **cd ~/neusomatic/test/**

[user@cn4192 test]$ **module load neusomatic**
[+] Loading neusomatic  0.2.1  on cn4192
[+] Loading singularity  on cn4192

[user@cn4192 test]$ **./run\_test.sh**
INFO 2019-10-24 17:23:17,745 __main__             Namespace(dbsnp_to_filter=None, del_merge_min_af=0, del_min_af=0.05, ensemble_tsv=None, filter_duplicate=False, first_do_without_qual=False, good_ao=10, ins_merge_min_af=0, ins_min_af=0.05, long_read=False, matrix_base_pad=7, matrix_width=32, max_dp=100000, merge_r=0.5, min_ao=1, min_dp=5, min_ev_frac_per_col=0.06, min_mapq=10, mode='call', normal_bam='../normal.bam', num_threads=1, reference='Homo_sapiens.GRCh37.75.dna.chromosome.22.fa', region_bed='../region.bed', restart=False, scan_alignments_binary='/opt/neusomatic/neusomatic/bin/scan_alignments', scan_maf=0.05, scan_window_size=2000, snp_min_af=0.05, snp_min_ao=10.0, snp_min_bq=20.0, truth_vcf=None, tsv_batch_size=50000, tumor_bam='../tumor.bam', work='work_standalone')
INFO 2019-10-24 17:23:17,745 preprocess           ----------------------Preprocessing------------------------
INFO 2019-10-24 17:23:17,747 preprocess           Scan tumor bam (and extracting quality scores).
INFO 2019-10-24 17:23:17,748 process_split_region Scan bam.
INFO 2019-10-24 17:23:17,748 scan_alignments      -------------------Scan Alignment BAM----------------------
[...]
INFO 2019-10-24 17:23:22,886 __main__             use_cuda: True
INFO 2019-10-24 17:23:22,886 call_neusomatic      -----------------Call Somatic Mutations--------------------
INFO 2019-10-24 17:23:22,886 call_neusomatic      PyTorch Version: 1.1.0
INFO 2019-10-24 17:23:22,886 call_neusomatic      Torchvision Version: 0.3.0
INFO 2019-10-24 17:23:22,931 call_neusomatic      GPU calling!
[...]
INFO 2019-10-24 17:23:44,911 postprocess          Postprocessing is Done.
### NeuSomatic stand-alone: SUCCESS! ###
### NeuSomatic ensemble: SUCCESS! ###

[user@cn4192 test]$ **exit**
exit
salloc.exe: Relinquishing job allocation 40260664

[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. neusomatic.sh). For example:



```

#!/bin/bash
module load neusomatic
preprocess.py \
        --mode call \
        --reference Homo_sapiens.GRCh37.75.dna.chromosome.22.fa \
        --region_bed region.bed \
        --tumor_bam tumor.bam \
        --normal_bam normal.bam \
        --work work_standalone \
        --scan_maf 0.05 \
        --min_mapq 10 \
        --snp_min_af 0.05 \
        --snp_min_bq 20 \
        --snp_min_ao 10 \
        --ins_min_af 0.05 \
        --del_min_af 0.05 \
        --num_threads 1 \
        --scan_alignments_binary $NEUSOMATIC_BIN/scan_alignments

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --partition=gpu --gres=gpu:k80:1,lscratch:10 --mem=20g -c14 neusomatic.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile (e.g. neusomatic.swarm). For example:



```

preprocess.py --mode call --region_bed region1.bed --tumor_bam tumor.bam   --normal_bam ../normal.bam [...]
preprocess.py --mode call --region_bed region2.bed --tumor_bam tumor.bam   --normal_bam ../normal.bam [...]
preprocess.py --mode call --region_bed region3.bed --tumor_bam tumor.bam   --normal_bam ../normal.bam [...]

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f neusomatic.swarm -g 20 -t 14 --partition=gpu --gres=gpu:k80:1,lscratch:10 --module neusomatic
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g *#*  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t *#* Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module neusomatic Loads the neusomatic module for each subjob in the swarm 
 | |
 | |
 | |








