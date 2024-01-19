

document.querySelector('title').textContent = 'MashMap: a fast adaptive algorithm for computing whole-genome homology maps ';
**MashMap: a fast adaptive algorithm for computing whole-genome homology maps** 


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



MashMap is an approximate algorithm for computing local alignment boundaries between long
DNA sequences. Given a minimum alignment length and an identity threshold, it computes the
desired alignment boundaries and identity estimates using kmer-based statistics, 
and maintains sufficient probabilistic guarantees on the output sensitivity.



### Reference:


* Chirag Jain, Sergey Koren, Alexander Dilthey, Adam M. Phillippy, Srinivas Aluru   

 *A Fast Adaptive Algorithm for Computing Whole-Genome Homology Maps*   

[Bioinformatics](https://academic.oup.com/bioinformatics/article/34/17/i748/5093242),
2018, **34**(17): i748-i756; doi: 10.1093/bioinformatics/bty597


Documentation
* [MashMap GitHub page](https://github.com/marbl/MashMap)


Important Notes
* Module Name: MashMap (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set
	+ **MASHMAP\_BIN**       executable directory
	+ **MASHMAP\_SRC**       source code directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=8g --cpus-per-task=16**
[user@cn3200 ~]$**module load MashMap** 
[+] Loading GSL 2.5 for GCC 7.2.0 ... 
[+] Loading MashMap  2.0 
[user@cn3200 ~]$ **mashmap -h**
-----------------
Mashmap is an approximate long read or contig mapper based on Jaccard
similarity
-----------------
Example usage:
$ mashmap -r ref.fa -q seq.fq [OPTIONS]
$ mashmap --rl reference_files_list.txt -q seq.fq [OPTIONS]

Available options
-----------------
-h, --help
    Print this help page

-r , --ref 
 an input reference file (fasta/fastq)[.gz]

--refList , --rl 
 a file containing list of reference files, one per line

-q , --query 
 an input query file (fasta/fastq)[.gz]

--ql , --queryList 
 a file containing list of query files, one per line

...

-t , --threads 
 count of threads for parallel execution [default : 1]

-o , --output 
 output file name [default : mashmap.out]
...

```

Download testing data from the PacBio's release of Human 54x long-read coverage dataset   

(<http://datasets.pacb.com/2014/Human54x/fastq.html>):

```

[user@cn3200 ~]$ **wget http://datasets.pacb.com/2013/Human10x/READS/2530572/0001/Analysis\_Results/m130929\_024849\_42213\_c100518541910000001823079209281311\_s1\_p0.1.subreads.fastq -O test\_seq.fq**

--2019-01-10 12:08:07--  http://datasets.pacb.com/2013/Human10x/READS/2530572/0001/Analysis_Results/m130929_024849_42213_c100518541910000001823079209281311_s1_p0.1.subreads.fastq
Resolving dtn04-e0 (dtn04-e0)... 10.1.200.240
Connecting to dtn04-e0 (dtn04-e0)|10.1.200.240|:3128... connected.
Proxy request sent, awaiting response... 200 OK
Length: 361764079 (345M) [application/octet-stream]
Saving to: ‘test_seq.fq’

100%[===============================================================>] 361,764,079 96.8MB/s   in 3.6s

2019-01-10 12:08:11 (96.2 MB/s) - ‘test_seq.fq’ saved [361764079/361764079]

```

Specify a human reference GRCh38 sequence:

```

[user@cn3200 ~]$ **ln -s /fdb/igenomes/Homo\_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa human\_GRCh38.fa**

```

Run the software on these data with 16 threads:

```

[user@cn3200 ~]$ **mashmap -r human\_GRCh38.fa -q test\_seq.fq -t 16**
>>>>>>>>>>>>>>>>>>
Reference = [human_GRCh38.fa]
Query = [test_seq.fq]
Kmer size = 16
Window size = 111
Segment length = 5000 (read split allowed)
Alphabet = DNA
Percentage identity threshold = 85%
Mapping output file = mashmap.out
Filter mode = 1 (1 = map, 2 = one-to-one, 3 = none)
Execution threads  = 16
>>>>>>>>>>>>>>>>>>
INFO, skch::Sketch::build, minimizers picked from reference = 52513331
INFO, skch::Sketch::index, unique minimizers = 17607584
INFO, skch::Sketch::computeFreqHist, Frequency histogram of minimizers = (1, 9361370) ... (116827, 1)
INFO, skch::Sketch::computeFreqHist, With threshold 0.001%, ignore minimizers occurring >= 5534 times during lookup.
INFO, skch::main, Time spent computing the reference index: 58.0907 sec
INFO, skch::Map::mapQuery, [count of mapped reads, reads qualified for mapping, total input reads] = [11978, 13497, 25249]
INFO, skch::main, Time spent mapping the query : 21.6042 sec
INFO, skch::main, mapping results saved in : mashmap.out

```

The result will be stored in the file 'mashmap.out'.   
  

End the interactive session:

```

[user@cn3200 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```





