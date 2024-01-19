

document.querySelector('title').textContent = 'uLTRA: a tool for splice alignment of long transcriptomic reads to a genome';
**uLTRA: a tool for splice alignment of long transcriptomic reads to a genome**


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



uLTRA implements an alignment method for long RNA sequencing reads 
based on a novel two-pass collinear chaining algorithm. uLTRA is guided by a database of exon annotations,
but it can also be used as a wrapper around minimap2 to align reads outside annotated regions.



### References:


* Kristoffer Sahlin and Veli Makinen   

Accurate spliced alignment of long RNA sequencing reads    

[Bioinformatics](https://academic.oup.com/bioinformatics/article/37/24/4643/6327681)  **37**(24), Pages 4643–4651 (2021)


Documentation
* [uLTRA GitHub page](https://github.com/ksahlin/ultra)


Important Notes
* Module Name: ultra (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set
	+ **ULTRA\_HOME**  installation directory
	+ **ULTRA\_BIN**       executable directory
	+ **ULTRA\_SRC**       source code directory
	+ **ULTRA\_DATA**  sample data directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=4g**
[user@cn0911 ~]$**module load ultra** 
[+] Loading singularity  3.10.5  on cn0911
[+] Loading Loading ultra  0.1
[user@cn0911 ~]$**cp $ULTRA\_DATA/\* .** 
[user@cn0911 ~]$**uLTRA pipeline SIRV\_genes.fasta SIRV\_genes\_C\_170612a.gtf reads.fa outfolder**
creating outfolder
total_flanks2: 90
total_flank_size 72816
total_unique_segment_counter 24396
total_segments_bad 25056
bad 460
total parts size: 24856
total exons size: 74187
min_intron: 20
Number of ref seqs in gff: 71
Number of ref seqs in fasta: 7
Filtering reads aligned to unindexed regions with minimap2
Running minimap2...
minimap2 done.
Done filtering. Reads filtered:0
Time for minimap2:0.10440826416015625 seconds.
Processing reads for NAM finding
Completed processing poly-A tails
Time for processing reads:0.006848812103271484 seconds.
RUNNING NAMFINDER USING COMMAND:
namfinder -k 10 -s 10 -l 10 -u 11 -C 500 -L 1000 -t 3 -S outfolder/refs_sequences.fa outfolder/reads_tmp.fa.gz 2> outfolder/namfinder_stderr.1 | gzip -1 --stdout > outfolder/seeds.txt.gz
Time for namfinder to find seeds:0.48410677909851074 seconds.
Starting aligning reads.
Wrote 0 batches of reads.
file_IO: Reading records done. Tot read: 4
file_IO: Written records in producer process: 0
Process: [0, 0, 0, 0, 0, 0, 0, 0]
Number of instances solved with quadratic collinear chainer solution: 0
Number of instances solved with n*log n collinear chainer solution: 0
Process: [0, 0, 0, 0, 0, 0, 0, 0]
Number of instances solved with quadratic collinear chainer solution: 0
Number of instances solved with n*log n collinear chainer solution: 0
Process: [3.8929521968943623, 4, 0, 0, 0, 0, 0, 0]
Number of instances solved with quadratic collinear chainer solution: 4
Number of instances solved with n*log n collinear chainer solution: 0
Wrote 1 batches of reads.
file_IO: Remainig written records after consumer join: 4
Done joining processes.
Time to align reads:0.07692599296569824 seconds.
Total mm2's primary alignments replaced with uLTRA: 2
Total mm2's alignments unmapped but mapped with uLTRA: 0
2 primary alignments had better score with uLTRA.
0 primary alignments had equal score with alternative aligner.
2 primary alignments had slightly better score with alternative aligner (typically ends bonus giving better scoring in ends, which needs to be implemented in uLTRA).
0 primary alignments had significantly better score with alternative aligner.
0 reads were unmapped with ultra but not by alternative aligner.
0 reads were not attempted to be aligned with ultra (unindexed regions), instead alternative aligner was used.
Time for selecting final best alignments (selecting the best of mm2's vs uLTRA's alignments):0.009923696517944336 seconds.
FSM : 4, NO_SPLICE : 0, Insufficient_junction_coverage_unclassified : 0, ISM/NIC_known : 0, NIC_novel : 0, NNC : 0
total alignment coverage: 3.8929521968943623
Deleting temporary files...
removed: outfolder/minimap2.sam
removed: outfolder/minimap2_errors.1
removed: outfolder/uLTRA_batch_0.stderr
removed: outfolder/uLTRA_batch_1.stderr
removed: outfolder/uLTRA_batch_2.stderr
removed: outfolder/reads_after_genomic_filtering.fasta
removed: outfolder/indexed.sam
removed: outfolder/unindexed.sam
removed: outfolder/refs_sequences.fa
Done.
[user@cn0911 ~]$**exit**

```
 
End the interactive session:

```

[user@cn0911 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```





