

document.querySelector('title').textContent = 'deepvariant on Biowulf';
deepvariant on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Interactive job - deeptrio](#int-trio)
[Batch job](#sbatch) 
[Swarm of jobs](#swarm) 
 |


 DeepVariant is an analysis pipeline that uses a deep neural network to call
genetic variants from next-generation DNA sequencing data *by converting
pileups from bam files to images and feeding them to a DNN based model*. 


Starting with version 1.5.0 a variant module with the `-deeptrio` is provided.
This variant is targeted at analyzing trios and duos. An example is provieded below


### References:


* R. Poplin, D. Newburger, J. Dijamco, N. Nguyen, D. Loy, S. S. Gross, C. Y. McLean, M. A. DePristo.
*Creating a universal SNP and small indel variant caller with deep neural networks*.
 BioRxiv, Dec 2017 [doi.org/10.1101/092890](https://www.biorxiv.org/content/early/2017/12/16/092890)


Documentation
* deepvariant Main Site: [GitHub](https://github.com/google/deepvariant)


Important Notes
* Module Name: deepvariant (see [the modules page](/apps/modules.html) for more information)
* Test data can be found in `${DEEPVARIANT_TEST_DATA}` (or `${DEEPTRIO_TEST_DATA}` for
 the deeptrio variant
* Please use GPU nodes only for step 2 of the Deepvariant pipeline (calling variants). The other
 steps don't benefit from GPU acceleration.



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Note that for the purposes of this
demonstration a GPU node is used for all three steps of the pipeline. For real applications, please only use GPU nodes for the `call_variant` step and
make sure to parallelize the `make_examples` step. Sample session:



```

[user@biowulf]$ **sinteractive --gres=gpu:p100:1,lscratch:50 --cpus-per-task=6**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn4224 are ready for job

[user@cn4224 ~]$ **cd /lscratch/$SLURM\_JOB\_ID**
[user@cn4224 ~]$ **module load deepvariant/1.5.0**
[user@cn4224 ~]$ **cp -r ${DEEPVARIANT\_TEST\_DATA:-none} input**
[user@cn4224 ~]$ **tree input**
input
|-- [user   3.7M]  NA12878_S1.chr20.10_10p1mb.bam
|-- [user   5.3K]  NA12878_S1.chr20.10_10p1mb.bam.bai
|-- [user    264]  test_nist.b37_chr20_100kbp_at_10mb.bed
|-- [user   5.6K]  test_nist.b37_chr20_100kbp_at_10mb.vcf.gz
|-- [user    197]  test_nist.b37_chr20_100kbp_at_10mb.vcf.gz.tbi
|-- [user    61M]  ucsc.hg19.chr20.unittest.fasta
|-- [user     23]  ucsc.hg19.chr20.unittest.fasta.fai
|-- [user   593K]  ucsc.hg19.chr20.unittest.fasta.gz
|-- [user     23]  ucsc.hg19.chr20.unittest.fasta.gz.fai
`-- [user    15K]  ucsc.hg19.chr20.unittest.fasta.gz.gzi

[user@cn4224 ~]$ **mkdir output**
[user@cn4224 ~]$ **REF=input/ucsc.hg19.chr20.unittest.fasta**
[user@cn4224 ~]$ **BAM=input/NA12878\_S1.chr20.10\_10p1mb.bam**

```

Extract images of pileups from the input bam file with `make_examples`, which
is single threaded and consumes 1-2GB of memory. Note that starting with version 1.1.0
it is necessary to specify `--channels insert_size` to make\_examples for WGS and
WES models. For other models other channels are needed.



```

[user@cn4224 ~]$ **make\_examples --mode calling \
 --ref "${REF}" \
 --reads "${BAM}" \
 --regions "chr20:10,000,000-10,010,000" \
 --channels insert\_size \
 --examples output/examples.tfrecord.gz**

```

Please note, that `make_examples` must be called with `--norealign_reads
--vsc_min_fraction_indels 0.12` flag for PacBio long reads.


This step can be parallelized with gnu parallel on a single machine, in which case multiple output files
are generated. Let's clean up the serial results and re-produce them in parallel.



```

[user@cn4224 ~]$ **rm output/examples.tfrecord.gz** # discard the output from the non-parallel run above
[user@cn4224 ~]$ **module load parallel**
[user@cn4224 ~]$ **N\_SHARDS=3**
[user@cn4224 ~]$ **mkdir -p logs**
[user@cn4224 ~]$ **seq 0 $((N\_SHARDS-1)) \
 | parallel -j ${SLURM\_CPUS\_PER\_TASK} --eta --halt 2 --joblog "logs/log" --res "logs" \
 make\_examples --mode calling \
 --ref "${REF}" \
 --reads "${BAM}" \
 --regions "chr20:10,000,000-10,010,000" \
 --examples "output/examples.tfrecord@${N\_SHARDS}.gz" \
 --channels insert\_size \
 --task {}**
[user@cn4224 ~]$ **tree output**
output
|-- [user   180K]  examples.tfrecord-00000-of-00003.gz
|-- [user   151K]  examples.tfrecord-00000-of-00003.gz.run_info.pbtxt
|-- [user   167K]  examples.tfrecord-00001-of-00003.gz
|-- [user   151K]  examples.tfrecord-00001-of-00003.gz.run_info.pbtxt
|-- [user   174K]  examples.tfrecord-00002-of-00003.gz
`-- [user   151K]  examples.tfrecord-00002-of-00003.gz.run_info.pbtxt

```

Call variants on the sharded output. The module is stored at /opt/models/model.chpt in
the singularity image containing deepvariant. To parallelize this step, multiple separate
jobs would have to be run. Note that since version 0.5.1, there are two models - one
for WGS data, one for WES. They are **/opt/wes/model.ckpt** and
**/opt/wes/model.ckpt**. For backwards compatibility with version 0.4.1,
**/opt/models/model.ckpt** points to the WGS model as well.



```

[user@cn4224 ~]$ **CALL\_VARIANTS\_OUTPUT="output/call\_variants\_output.tfrecord.gz"**
[user@cn4224 ~]$ **#< 0.9.0 MODEL="/opt/wgs/model.ckpt"**
[user@cn4224 ~]$ **MODEL="/opt/models/wgs/model.ckpt"**
[user@cn4224 ~]$ **call\_variants \
 --outfile "${CALL\_VARIANTS\_OUTPUT}" \
 --examples "output/examples.tfrecord@${N\_SHARDS}.gz" \
 --checkpoint "${MODEL}"**

```

Finally, convert the output to VCF format



```

[user@cn4224 ~]$ **postprocess\_variants \
 --ref "${REF}" \
 --infile "${CALL\_VARIANTS\_OUTPUT}" \
 --outfile "output/example.vcf.gz"**
[user@cn4224 ~]$ **tree output**
output
|-- [user   4.1K]  call_variants_output.tfrecord.gz
|-- [user   180K]  examples.tfrecord-00000-of-00003.gz
|-- [user   151K]  examples.tfrecord-00000-of-00003.gz.run_info.pbtxt
|-- [user   166K]  examples.tfrecord-00001-of-00003.gz
|-- [user   151K]  examples.tfrecord-00001-of-00003.gz.run_info.pbtxt
|-- [user   174K]  examples.tfrecord-00002-of-00003.gz
|-- [user   151K]  examples.tfrecord-00002-of-00003.gz.run_info.pbtxt
`-- [user   2.2K]  example.vcf.gz

[user@cn4224 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```

With version 0.5.1, deepvariant gained the ability to generate gVCF output.


In version 0.9.0 deepvariant added a single command (`run_deepvariant`) to
run the whole pipeline. However, since step 1 and 3 can be run without GPU it is more efficient
to run the steps separately


Interactive job - deeptrio
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Deeptrio is build on top of deepvariant and follows the same steps. Deeptrio considers genetic inheritance 
when calling variants and produces joint examples from all samples to ensure that variant calls are made
for the same sites in all samples. It is available as a variant version starting with 1.5.0.



```

[user@biowulf]$ **sinteractive --gres=gpu:p100:1,lscratch:50 --cpus-per-task=6**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn4224 are ready for job

[user@cn4224 ~]$ **cd /lscratch/$SLURM\_JOB\_ID**
[user@cn4224 ~]$ **module load deepvariant/1.5.0-deeptrio**
[user@cn4224 ~]$ **cp -r ${DEEPTRIO\_TEST\_DATA:-none} input**
[user@cn4224 ~]$ **tree input**
input
|-- [user   2.9G]  GRCh38_no_alt_analysis_set.fasta
|-- [user   7.6K]  GRCh38_no_alt_analysis_set.fasta.fai
|-- [user   1.3M]  HG002.chr20.10_10p1mb.bam
|-- [user   6.7K]  HG002.chr20.10_10p1mb.bam.bai
|-- [user    11M]  HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
|-- [user   149M]  HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
|-- [user   1.6M]  HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi
|-- [user   1.4M]  HG003.chr20.10_10p1mb.bam
|-- [user   6.7K]  HG003.chr20.10_10p1mb.bam.bai
|-- [user    13M]  HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
|-- [user   140M]  HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
|-- [user   1.6M]  HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi
|-- [user   1.3M]  HG004.chr20.10_10p1mb.bam
|-- [user   6.7K]  HG004.chr20.10_10p1mb.bam.bai
|-- [user    12M]  HG004_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
|-- [user   142M]  HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
`-- [user   1.6M]  HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi

[user@cn4224 ~]$ **mkdir output**

```


Make examples

```

[user@cn4224 ~]$ **module load parallel**
[user@cn4224 ~]$ **N\_SHARDS=$SLURM\_CPUS\_PER\_TASK**
[user@cn4224 ~]$ **mkdir -p logs output**
[user@cn4224 ~]$ **seq 0 $((N\_SHARDS-1)) \
 | parallel -j ${SLURM\_CPUS\_PER\_TASK} --eta --halt 2 --joblog "logs/log" --res "logs" \
 make\_examples --mode calling \
 --ref input/GRCh38\_no\_alt\_analysis\_set.fasta \
 --reads\_parent1 "input/HG003.chr20.10\_10p1mb.bam" \
 --sample\_name\_parent1 "HG003" \
 --reads\_parent2 "input/HG004.chr20.10\_10p1mb.bam" \
 --sample\_name\_parent2 "HG004" \
 --reads "input/HG002.chr20.10\_10p1mb.bam" \
 --sample\_name "HG002" \
 --regions "chr20:10,000,000-10,010,000" \
 --examples "output/intermediate\_results\_dir/make\_examples.tfrecord@${N\_SHARDS}.gz" \
 --channels insert\_size \
 --gvcf "output/intermediate\_results\_dir/gvcf.tfrecord@6.gz" \
 --pileup\_image\_height\_child "60" \
 --pileup\_image\_height\_parent "40" \
 --task {}**

```

Call variants on the sharded output. The module is stored at /opt/models/model.chpt in
the singularity image containing deepvariant. To parallelize this step, multiple separate
jobs would have to be run - for example a job for each member of the trio



```

[user@cn4224 ~]$ **for n in parent1 parent2 child ; do
 call\_variants \
 --outfile output/intermediate\_results\_dir/call\_variants\_output\_${n}.tfrecord.gz \
 --examples output/intermediate\_results\_dir/make\_examples\_${n}.tfrecord@${N\_SHARDS}.gz \
 --checkpoint "/opt/models/deeptrio/wgs/parent/model.ckpt"
 done**

```

Finally, convert the output to VCF format



```

[user@cn4224 ~]$ **declare -Asamples=( [parent1]=HG003 [parent2]=HG004 [child]=HG002 )**
[user@cn4224 ~]$ **for n in parent1 parent2 child ; do
 postprocess\_variants \
 --ref input/GRCh38\_no\_alt\_analysis\_set.fasta \
 --infile output/intermediate\_results\_dir/call\_variants\_output\_${n}.tfrecord.gz \
 --outfile output/${samples[$n]}.output.vcf.gz \
 --nonvariant\_site\_tfrecord\_path output/intermediate\_results\_dir/gvcf\_${n}.tfrecord@${N\_SHARDS}.gz \
 --gvcf\_outfile output/${samples[$n]}.g.vcf.gz
 done**
[user@cn4224 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file for each step (e.g. deepvariant\_step[123].sh) similar to the following. Note
that for simplicity there is little error checking in the example scripts. Note also that only step2
makes use of a GPU and that the current version of DeepVariant cannot scale to more than 1 GPU.



```

#! /bin/bash
# this is deepvariant_step1.sh

module load deepvariant/1.4.0 || exit 1
module load parallel

wd=$PWD

[[ -d /lscratch/${SLURM_JOB_ID} ]] && cd /lscratch/${SLURM_JOB_ID} || exit 1
rm -rf input output
mkdir input output $wd/logs-parallel-$SLURM_JOB_ID
cp "${DEEPVARIANT_TEST_DATA}"/* input || exit 1

REF=input/ucsc.hg19.chr20.unittest.fasta
BAM=input/NA12878_S1.chr20.10_10p1mb.bam
N_SHARDS=12

seq 0 $((N_SHARDS-1)) \
    | parallel -P ${SLURM_CPUS_PER_TASK} --halt 2 \
        --joblog "$wd/logs-parallel-$SLURM_JOB_ID/log" --res "$wd/logs-parallel-$SLURM_JOB_ID" \
      make_examples --mode calling \
        --ref "${REF}" \
        --reads "${BAM}" \
        --regions "chr20:10,000,000-10,010,000" \
        --examples output/examples.tfrecord@${N_SHARDS}.gz\
        --channels insert_size \
        --task {} \
|| exit 1

cp -r output $wd

```


```

#! /bin/bash
# this is deepvariant_step2.sh

module load deepvariant/1.4.0 || exit 1

#< 0.9.0 MODEL="/opt/wgs/model.ckpt"
MODEL="/opt/models/wgs/model.ckpt"
N_SHARDS=12
CALL_VARIANTS_OUTPUT="output/call_variants_output.tfrecord.gz"

call_variants \
 --outfile "${CALL_VARIANTS_OUTPUT}" \
 --examples "output/examples.tfrecord@${N_SHARDS}.gz" \
 --checkpoint "${MODEL}" \
|| exit 1

```


```

#! /bin/bash
# this is deepvariant_step3.sh

module load deepvariant/1.4.0 || exit 1

REF=input/ucsc.hg19.chr20.unittest.fasta
CALL_VARIANTS_OUTPUT="output/call_variants_output.tfrecord.gz"

postprocess_variants \
  --ref "${REF}" \
  --infile "${CALL_VARIANTS_OUTPUT}" \
  --outfile "output/example.vcf.gz" \
|| exit 1

```

Submit these jobs using the Slurm [sbatch](/docs/userguide.html) command. In
this example, job dependencies are used to tie the jobs together.



```

biowulf$ **sbatch --cpus-per-task=6 --mem-per-cpu=2g --gres=lscratch:10 deepvariant\_step1.sh**
1111
biowulf$ **sbatch --dependency=afterany:1111 --cpus-per-task=6 \
 --gres=lscratch:10,gpu:k80:1 --partition=gpu deepvariant\_step2.sh**
1112
biowulf$ **sbatch --dependency=afterany:1112 deepvariant\_step3.sh**
1113

```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resources.
Create a swarmfile for the first step of the pipeline (e.g. deepvariant.swarm). For example:



```

make_examples --mode calling --ref /path/to/genome.fa --reads "sample1.bam" \
    --regions "chr1" --examples sample1_chr1.tfrecord.gz --channels insert_size
make_examples --mode calling --ref /path/to/genome.fa --reads "sample1.bam" \
    --regions "chr2" --examples sample1_chr2.tfrecord.gz --channels insert_size
make_examples --mode calling --ref /path/to/genome.fa --reads "sample1.bam" \
    --regions "chr3" --examples sample1_chr3.tfrecord.gz --channels insert_size

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f deepvariant.swarm [-g #] --module deepvariant
```

where


|  |  |  |  |
| --- | --- | --- | --- |
| -g #  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | --module deepvariant  Loads the deepvariant module for each subjob in the swarm
 | |
 | |








