

document.querySelector('title').textContent = 'Liger2LiGer: Nanopore chimera splitting/detection';
**Liger2LiGer: Nanopore chimera splitting/detection**


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



Liger2LiGer is a Nanopore chimera splitting/detection tool. 
Automated end-to-end chimera detection and dataset evaluation starting from fastq.



Documentation
* [Liger2LiGer GitHub page](https://github.com/rlorigro/Liger2LiGer)


Important Notes
* Module Name: liger2liger (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set
	+ **L2L\_HOME**  installation directory
	+ **L2L\_BIN**       executable directory
	+ **L2L\_SRC**       source code directory
	+ **L2L\_DATA**  sample data directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=4g -c4**
[user@cn0911 ~]$**module load liger2liger** 
[+] Loading singularity  3.10.5  on cn0911
[+] Loading liger2liger  20230901
[user@cn0911 ~]$**evaluate\_chimeras.py -h**
usage: evaluate_chimeras.py [-h] --ref REF --fastq FASTQ [--dry] [--n_threads N_THREADS]

optional arguments:
  -h, --help            show this help message and exit
  --ref REF, -r REF     Path of reference fasta to use for alignment
  --fastq FASTQ, -i FASTQ
                        comma-separated list of fastq file paths to be aligned
  --dry                 echo the commands instead of executing them
  --n_threads N_THREADS
                        How many threads to use for minimap2 alignment

[user@cn0911 ~]$**extract\_reads\_from\_fastq.py -h**
usage: extract_reads_from_fastq.py [-h] --fastq FASTQ --ids IDS

optional arguments:
  -h, --help     show this help message and exit
  --fastq FASTQ  path of file containing FASTA/FASTQ sequence
  --ids IDS      path of file containing 1 id per line to be queried

[user@cn0911 ~]$**filter\_paf\_by\_read\_name.py -h**
usage: filter_paf_by_read_name.py [-h] --paf PAF --names NAMES

optional arguments:
  -h, --help     show this help message and exit
  --paf PAF      path of PAF file to be filtered
  --names NAMES  path of file containing a list of read ids, one on each line
[user@cn0911 ~]$**generate\_chimer\_stats.py**
usage: generate_chimer_stats.py [-h] --input INPUT [--output_dir OUTPUT_DIR] [--dry]

optional arguments:
  -h, --help            show this help message and exit
  --input INPUT, -i INPUT
                        Comma separated paths of length txt files containing lengths of chimers and non-chimers (one length [in
                        bp] per line). Can run any number of pairs of txt files as long as they have the suffix
                        non_chimer_lengths.txt or chimer_lengths.txt and pairs share a prefix.
  --output_dir OUTPUT_DIR, -o OUTPUT_DIR
                        Optionally create a directory to store the output.
  --dry, -d             Don't create any files, report output paths only.

```

Run liger2liger on sample data:

```

[user@cn0911 ~]$**mkdir /data/$USER/L2L && cd /data/$USER/L2L**
[user@cn0911 ~]$**cp $L2L\_DATA/\* .** 
[user@cn0911 ~]$**ln -s /fdb/igenomes/Drosophila\_melanogaster/UCSC/dm6/Sequence/Chromosomes/chrX.fa**
[user@cn0911 ~]$**evaluate\_chimeras.py --ref chrX.fa --fastq test.fastq --n\_threads 2** 
STARTING: test_VS_chrX
Found executable: /opt/conda/envs/l2l/Liger2LiGer/build/filter_chimeras_from_alignment
minimap2 -x map-ont --secondary=no -n 10 -K 10g -k 17 -t 2 chrX.fa test.fastq
[M::mm_idx_gen::0.863*1.00] collected minimizers
[M::mm_idx_gen::1.122*1.23] sorted minimizers
[M::main::1.122*1.23] loaded/built the index for 1 target sequence(s)
[M::mm_mapopt_update::1.197*1.22] mid_occ = 30
[M::mm_idx_stat] kmer size: 17; skip: 10; is_hpc: 0; #seq: 1
[M::mm_idx_stat::1.255*1.21] distinct minimizers: 3830080 (96.13% are singletons); average occurrences: 1.118; average spacing: 5.500
[M::worker_pipeline::1.276*1.21] mapped 350 sequences
[M::main] Version: 2.17-r941
[M::main] CMD: minimap2 -x map-ont --secondary=no -n 10 -K 10g -k 17 -t 2 chrX.fa test.fastq
[M::main] Real time: 1.286 sec; CPU: 1.551 sec; Peak RSS: 0.216 GB
Writing chimeric reads to file: "/vf/users/user/L2L/test_VS_chrX/test_VS_chrX.chimeric_reads.txt"
Writing non-chimeric reads to file: "/vf/users/user/L2L/test_VS_chrX/test_VS_chrX.non_chimeric_reads.txt"
Writing chimeric lengths to file: "/vf/users/user/L2L/test_VS_chrX/test_VS_chrX.chimer_lengths.txt"
Writing non-chimeric lengths to file: "/vf/users/user/L2L/test_VS_chrX/test_VS_chrX.non_chimer_lengths.txt"
/vf/users/user/L2L/test_VS_chrX/test_VS_chrX.non_chimer_lengths.txt
test_VS_chrX True True
/vf/users/user/L2L/test_VS_chrX/test_VS_chrX.chimer_lengths.txt
test_VS_chrX False True
test
N50     2625.0
Writing outputs to: /vf/users/user/L2L/test_VS_chrX/results_09_27_2023_09:03:52.csv
Writing outputs to: results_09_27_2023_09:03:52.csv
[user@cn0911 ~]$**exit**

```
 
End the interactive session:

```

[user@cn0911 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```





