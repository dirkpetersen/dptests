

document.querySelector('title').textContent = 'mbg on Biowulf';
mbg on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
 |



From the repository:




> 
> Minimizer based sparse de Bruijn Graph constructor. Homopolymer
> compress input sequences, pick syncmers from hpc-compressed sequences,
> connect syncmers with an edge if they are adjacent in a read, unitigify and
> homopolymer decompress. Suggested input is PacBio HiFi/CCS reads, or ONT duplex
> reads. May or may not work with Illumina reads. Not suggested for PacBio CLR or
> regular ONT reads
> 


### References:


* M. Rautiainen, T. Marschall. *MBG: Minimizer-based Sparse de Bruijn Graph Construction*.
 Genome Biology (2020). [PubMed](https://pubmed.ncbi.nlm.nih.gov/33475133/) | 
 [PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/pmid/33475133/) | 
 [Journal](https://academic.oup.com/bioinformatics/article/37/16/2476/6104877?login=true)


Documentation
* mbg on [GitHub](https://github.com/maickrau/MBG)


Important Notes
* Module Name: mbg (see [the modules page](/apps/modules.html) for more information)
* MBG is a multithreaded tool. Please match the number of allocated CPUs to the number of threads
* Example files in `$MBG_TEST_DATA`



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --gres=lscratch:10 --cpus-per-task=2 --mem=3g**
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144]$ **cd /lscratch/$SLURM\_JOB\_ID**
[user@cn3144]$ **module load mbg**
[user@cn3144]$ **cp ${MBG\_TEST\_DATA:-none}/SRR10971019.fasta .**
[user@cn3144]$ **MBG -t $SLURM\_CPUS\_PER\_TASK -i SRR10971019.fasta -o SRR10971019\_graph.gfa -k 1501 -w 1450 -a 1 -u 3**
Parameters: k=1501,w=1450,a=1,u=3,t=2,r=0,R=0,hpcvariantcov=0,errormasking=hpc,endkmers=no,blunt=no,keepgaps=no,guesswork=no,cache=no
Collecting selected k-mers
Reading sequences from SRR10971019.fasta
1210730 total selected k-mers in reads
265228 distinct selected k-mers in reads
Unitigifying
Filtering by unitig coverage
3513 distinct selected k-mers in unitigs after filtering
Getting read paths
Reading sequences from SRR10971019.fasta
Building unitig sequences
Reading sequences from SRR10971019.fasta
Writing graph to SRR10971019_graph.gfa
selecting k-mers and building graph topology took 19,594 s
unitigifying took 0,81 s
filtering unitigs took 0,4 s
getting read paths took 19,186 s
building unitig sequences took 36,835 s
forcing edge consistency took 0,24 s
writing the graph and calculating stats took 0,94 s
nodes: 567
edges: 730
assembly size 5346906 bp, N50 29122
approximate number of k-mers ~ 4495839

[user@cn3144]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. mbg.sh), which uses the input file 'mbg.in'. For example:



```

#!/bin/bash
module load mbg/1.0.11
cp ${MBG_TEST_DATA:-none}/SRR10971019.fasta .
MBG -t $SLURM_CPUS_PER_TASK -i SRR10971019.fasta -o SRR10971019_graph.gfa -k 1501 -w 1450 -a 1 -u 3

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch [--cpus-per-task=#] [--mem=#] mbg.sh
```







