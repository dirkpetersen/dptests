

document.querySelector('title').textContent = 'deepsignal on Biowulf';
deepsignal on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
 |



From the deepsignal home page



> A deep-learning method for detecting DNA methylation state from Oxford
>  Nanopore sequencing reads.


### References:


* P. Ni, N. Huang, Z. Zhang, D. Wang, F. Liang, Y. Miao, C. Xiao, F. Luo, and J. Wang. 
 *DeepSignal: detecting DNA methylation state from Nanopore sequencing reads using 
 deep-learning.*, Bioinformatics 2019, 35:4586-4595.
 [PubMed](https://pubmed.ncbi.nlm.nih.gov/30994904/) | 
 [Journal](https://academic.oup.com/bioinformatics/article/35/22/4586/5474907)



Documentation
* deepsignal [on GitHub](https://github.com/bioinfomaticsCSU/deepsignal)


Important Notes
* Module Name: deepsignal (see [the modules page](/apps/modules.html) for more information)
* guppy used in the example and deepsignal make use of GPUs
* Example files in `$DEEPSIGNAL_TEST_DATA`



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=10g -c6 --gres=lscratch:50,gpu:p100:1**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144]$ **cd /lscratch/$SLURM\_JOB\_ID**
[user@cn3144]$ **module load guppy deepsignal/2-0.1.3**
[+] Loading singularity  3.7.1  on cn2384
[+] Loading guppy 4.2.2  ...
[+] Loading deepsignal  0.1.8

[user@cn3144]$ **tar -xzf ${DEEPSIGNAL\_TEST\_DATA:-none}/fast5s.sample.tar.gz**
[user@cn3144]$ **ll**
total 13M
drwxr-xr-x 3 user group 548K Dec 17  2018 fast5s.al
-rw-r--r-- 1 user group  12M Dec 17  2018 GCF_000146045.2_R64_genomic.fna

```

To call modifications, the raw fast5 files should be basecalled (Guppy or
Albacore) and then be re-squiggled by tombo. At last, modifications of
specified motifs can be called by deepsignal. The following are commands to
call 5mC in CG contexts from the example data:



```

[user@cn3144]$ **guppy\_basecaller -x cuda:all \
 -i fast5s.al -r -s fast5s.al.guppy --config dna\_r9.4.1\_450bps\_hac\_prom.cfg**
ONT Guppy basecalling software version 4.2.2+effbaf8
config file:        /opt/ont/guppy/data/dna_r9.4.1_450bps_hac_prom.cfg
model file:         /opt/ont/guppy/data/template_r9.4.1_450bps_hac_prom.jsn
input path:         fast5s.al
save path:          fast5s.al.guppy
chunk size:         2000
chunks per runner:  1248
records per file:   4000
num basecallers:    4
gpu device:         cuda:all
kernel path:
runners per device: 12

Found 3621 fast5 files to process.
Init time: 6689 ms

0%   10   20   30   40   50   60   70   80   90   100%
|----|----|----|----|----|----|----|----|----|----|
***************************************************
Caller time: 23257 ms, Samples called: 53658471, samples/s: 2.3072e+06
Finishing up any open output files.
Basecalling completed successfully.

[user@cn3144]$ **tombo preprocess annotate\_raw\_with\_fastqs \
 --fast5-basedir fast5s.al \
 --fastq-filenames fast5s.al.guppy/{pass,fail}/\*.fastq \
 --overwrite \
 --processes 6 \
 --sequencing-summary-filenames fast5s.al.guppy/sequencing\_summary.txt**
[14:45:09] Getting read filenames.
[14:45:09] Parsing sequencing summary files.
[14:45:09] Annotating FAST5s with sequence from FASTQs.
100%|███████████████████████████████████████████████████████████████| 3621/3621 
[14:45:11] Added sequences to a total of 3621 reads.

[user@cn3144]$ **tombo resquiggle fast5s.al \
 GCF\_000146045.2\_R64\_genomic.fna \
 --processes 6 \
 --corrected-group RawGenomeCorrected\_001 \
 --overwrite**
[14:47:15] Loading minimap2 reference.
[14:47:16] Getting file list.
[14:47:16] Loading default canonical ***** DNA ***** model.
[14:47:16] Re-squiggling reads (raw signal to genomic sequence alignment).
100%|█████████████████████████████████████████████████████████████████| 3621/3621
[14:49:08] Final unsuccessful reads summary (4.7% reads unsuccessfully processed; 169 total reads):
     3.2% (115 reads) : Alignment not produced
     1.2% ( 44 reads) : Poor raw to expected signal matching (revert with `tombo filter clear_filters`)
     0.3% ( 10 reads) : Read event to sequence alignment extends beyond bandwidth

[14:49:08] Saving Tombo reads index to file.

```

Notes:


* Some versions of deepsignal on biowulf include their own tombo. Others use the
 tombo module.
* Tombo does not take advantage of GPU. In the example above tombo is run in the same
 session as deepsignal only for convenience. For real data please run tombo
 as a separate job without GPU for efficiency.


The available models can be checked with `ls /usr/local/apps/deepsignal/VER/model` or `/usr/local/apps/deepsignal/VER/share/model`. Because deepsignal is installed as a container the path used for deepsignal is `/model/NAME`.



```

[user@cn3144]$ **deepsignal2 call\_mods \
 --nproc $SLURM\_CPUS\_PER\_TASK \
 --input\_path fast5s.al/ \
 --model\_path /model/model.dp2.CG.R9.4\_1D.human\_hx1.bn17\_sn16.both\_bilstm.b17\_s16\_epoch4.ckpt \
 --result\_file fast5s.al.CpG.call\_mods.tsv \
 --reference\_path GCF\_000146045.2\_R64\_genomic.fna \
 --corrected\_group RawGenomeCorrected\_001**
[...snip...]

[user@cn3144]$ **call\_modification\_frequency.py \
 --input\_path fast5s.al.CpG.call\_mods.tsv \
 --result\_file fast5s.al.CpG.call\_mods.frequency.tsv**
get 1 input file(s)..
reading the input files..
162916 of 162916 samples used..
writing the result..

[user@cn3144]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. deepsignal.sh). In this example we assume
that the base calling and resquiggle has already been completed and we just run
the call\_mods step



```

#!/bin/bash
module load deepsignal/2-0.1.3
deepsignal2 call_mods \
    --nproc $SLURM_CPUS_PER_TASK \
    --input_path fast5s.al/ \
    --model_path /model/model.dp2.CG.R9.4_1D.human_hx1.bn17_sn16.both_bilstm.b17_s16_epoch4.ckpt \
    --result_file fast5s.al.CpG.call_mods.tsv \
    --reference_path GCF_000146045.2_R64_genomic.fna \
    --corrected_group RawGenomeCorrected_001

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=6 --mem=20g --gres=gpu:p100:1,lscratch:10 deepsignal.sh
```







