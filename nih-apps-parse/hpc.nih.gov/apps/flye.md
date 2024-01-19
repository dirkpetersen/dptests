

document.querySelector('title').textContent = 'flye on Biowulf';
flye on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
 |



> 
> Flye is a de novo assembler for long and noisy reads, such as those produced by
> PacBio and Oxford Nanopore Technologies. The algorithm uses an A-Bruijn graph
> to find the overlaps between reads and does not require them to be
> error-corrected. After the initial assembly, Flye performs an extra repeat
> classification and analysis step to improve the structural accuracy of the
> resulting sequence. The package also includes a polisher module, which produces
> the final assembly of high nucleotide-level quality.
> 


Flye replaces abruijn and does provide a `abruijn` script for
backwards compatibility.



A 5Mb bacterial genome with ~80x coverage was assembled on
one of our compute nodes (6GB memory; 16 CPUs) in about 30min. A ~150 Mb D.
melanogaster genome was assembled in 13h (100GB memory; 32 CPUs).



### References:


* Yu Lina, Jeffrey Yuana, Mikhail Kolmogorova, Max W. Shena, Mark Chaissonb, and Pavel A. Pevzner.
 *Assembly of long error-prone reads using de Bruijn graphs* PNAS 2016, 27:E8396-E8405
 [PubMed](https://www.ncbi.nlm.nih.gov/pubmed/27956617) | 
 PMC  | 
 [Journal](http://www.pnas.org/content/113/52/E8396)


Documentation
* [GitHub](https://github.com/fenderglass/ABruijn)


Important Notes
* Module Name: flye (see [the modules page](/apps/modules.html) for more information)
* flye is a multithreaded application
* Example files can be found in `$FLYE_TEST_DATA`



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:


Note that `--genome-size` is optional since version 2.8



```


[user@biowulf]$ **sinteractive --mem=14g --cpus-per-task=16 --gres=lscratch:10**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job
[user@cn3144 ~]$ **module load flye**
[user@cn3144 ~]$ **cd /lscratch/$SLURM\_JOB\_ID**
[user@cn3144 ~]$ **zcat $FLYE\_TEST\_DATA/SRR1284073\_gt10k.fasta.gz > SRR1284073\_gt10k.fasta**
[user@cn3144 ~]$ **flye -t $SLURM\_CPUS\_PER\_TASK --pacbio-raw SRR1284073\_gt10k.fasta \
 -o assembly\_ecoli --genome-size 5m**
[2018-04-03 16:08:46] INFO: Running Flye 2.3.3-g0fc9012
[2018-04-03 16:08:46] INFO: Assembling reads
[2018-04-03 16:08:46] INFO: Running with k-mer size: 15
[2018-04-03 16:08:46] INFO: Reading sequences
[2018-04-03 16:08:53] INFO: Reads N50/90: 17480 / 11579
[2018-04-03 16:08:53] INFO: Selected minimum overlap 5000
[2018-04-03 16:08:53] INFO: Expected read coverage: 80
[...snip...]
[2018-04-03 16:31:44] INFO: Correcting bubbles
0% 10% 20% 30% 40% 50% 60% 70% 80% 90% 100%
[2018-04-03 16:36:05] INFO: Assembly statistics:

        Total length:   4636855
        Contigs:        1
        Scaffolds:      1
        Scaffolds N50:  4636855
        Largest scf:    4636855
        Mean coverage:  52

[2018-04-03 16:36:05] INFO: Final assembly: /lscratch/46116226/assembly_ecoli/scaffolds.fasta

[user@cn3144 ~]$ **ll assembly\_ecoli**
total 9.1M
drwxrwxr-x 2 user group 4.0K Apr  3 12:13 0-assembly
drwxrwxr-x 2 user group 4.0K Apr  3 12:21 1-consensus
drwxrwxr-x 2 user group 4.0K Apr  3 12:23 2-repeat
drwxrwxr-x 2 user group 4.0K Apr  3 12:36 3-polishing
-rw-rw-r-- 1 user group  193 Apr  3 12:23 assembly_graph.dot
-rw-rw-r-- 1 user group   79 Apr  3 12:36 assembly_info.txt
-rw-rw-r-- 1 user group 4.5M Apr  3 12:36 contigs.fasta
-rw-rw-r-- 1 user group  21K Apr  3 12:36 flye.log
-rw-rw-r-- 1 user group   26 Apr  3 12:36 flye.save
-rw-rw-r-- 1 user group 4.5M Apr  3 12:36 scaffolds.fasta


[user@cn3144 ~]$ # copy back to data
[user@cn3144 ~]$ **cp -r assembly\_ecoli /data/$USER/badbadproject**
[user@cn3144 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. flye.sh) similar to the following:



```

#! /bin/bash
ml flye || exit 1

cd /lscratch/$SLURM_JOB_ID
cp /data/users/some/where/reads.fa.gz .
flye -t $SLURM_CPUS_PER_TASK --pacbio-raw reads.fasta.gz -o assembly_dmelanogaster --genome-size 150m
mv assembly_dmelanogaster /data/$USER/badbadproject

```

This particular example made use of data set
[SRX499318](https://www.ncbi.nlm.nih.gov/sra/SRX499318[accn]) filtered
to reads >14k length resulting in a 90x coverage of the ~150Mb *D. melanogaster* genome.


Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```

sbatch --mem=120g --cpus-per-task=32 --gres=lscratch:300 flye.batch --time=1-00:00:00

```

This job ran for ~13h and used up to 100GB of memory. Here is the
profile of memory and running threads for this assembly:



![resource usage profile](/images/flye_dmelanogaster.png)

The final result was an assembly of 137Mb with 357 contigs and a scaffold N50 of 6.34Mb.








