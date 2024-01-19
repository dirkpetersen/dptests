

document.querySelector('title').textContent = ' HiNT on Biowulf';
 HiNT on Biowulf


|  |
| --- |
| 
Quick Links
[Interactive Jobs](#int)
[Batch job on Biowulf](#batch)
[Using Swarm](#swarm)
[Documentation](#doc)
 |


HiNT (Hi-C for copy Number variation and Translocation detection), a computational method to detect CNVs and Translocations from Hi-C data. HiNT has three main components: HiNT-PRE, HiNT-CNV, and HiNT-TL. HiNT-PRE preprocesses Hi-C data and computes the contact matrix, which stores contact frequencies between any two genomic loci; both HiNT-CNV and HiNT-TL starts with HI-C contact matrix, predicts copy number segments, and inter-chromosomal translocations, respectively



Documentation
<https://github.com/parklab/HiNT>


Important Notes
* Module Name: hint (see [the modules page](/apps/modules.html) for more information)
* see /fdb/hint for reference, index, matrices, and example files
* Programs are multithreaded.


Submitting an interactive job

Allocate an interactive session and run the interactive job there.



```
[biowulf]$ **sinteractive --mem=40g --cpus-per-task=16**
salloc.exe: Granted job allocation 789523
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn0135 are ready for job

[cn0135]$ **cd /data/$USER/**

[cn0135]$ **module load hint**

[cn0135]$ **cp /fdb/hint/testData/ .** 

[cn0135]$ **cd testData**

[cn0135]$ **hint pre -d $TESTDIR/test.bam \
--refdir $REFDIR/hg19 \
--informat bam \
--outformat juicer \
-g hg19 \
-n test \
-o /data/$USER/out \
--samtoolspath $SAMTOOLS \
-p $SLURM\_CPUS\_PER\_TASK \
--pairtoolspath $PAIRTOOLS \
--juicerpath $JUICERPATH/scripts/juicer\_tools\_1.22.01.jar**
[14:42:34] Argument List:
[14:42:34] Hi-C data = /, f, d, b, /, h, i, n, t, /, t, e, s, t, D, a, t, a, /, t, e, s, t, ., b, a, m
[14:42:34] Input format = bam
[14:42:34] Output format = juicer
[14:42:34] Genome = hg19
...
...
Calculating norms for zoom BP_10000
Calculating norms for zoom BP_5000
Writing expected
Writing norms
Finished writing norms

or

[cn0135]$ **hint pre \
-d $TESTDIR/TestSub\_1.fq.gz,$TESTDIR/TestSub\_2.fq.gz \
-a $BWA \
-i $BWAINDEXDIR/hg19/hg19.fa \
--refdir $REFDIR/hg19 \
--informat fastq \
--outformat cooler \
-g hg19 \
-n test \
-o /data/$USER/out \
--samtoolspath $SAMTOOLS \
-p $SLURM\_CPUS\_PER\_TASK \
--pairtoolspath $PAIRTOOLS \
--coolerpath $COOLER**

or

[cn0135]$ **hint tl \
-m test.hic \
-f juicer \
--refdir $REFDIR/hg19 \
--backdir $MATRICESDIR/hg19 \
-g hg19 \
-n test \
-c 0.05 \
--ppath $PAIRIX \
-p $SLURM\_CPUS\_PER\_TASK \
-o testout**

[cn0135]$ **exit**
salloc.exe: Job allocation 789523 has been revoked.
[biowulf]$
```


Note: this job allocates 10 GB of memory and automatically assign the number of cpus allocated to the variable $SLURM\_CPUS\_PER\_TASK. 

The test takes less than 30 minutes.



Submitting a single batch job
1. Create a script file (myscript) similar to the one below.

 

```
#!/bin/bash 

cd /data/$USER/testData
module load hint
hint pre -d $TESTDIR/test.bam \
--refdir $REFDIR/hg19 \
--informat bam \
--outformat juicer \
-g hg19 \
-n test \
-o /data/$USER/out \
--samtoolspath $SAMTOOLS \
-p $SLURM_CPUS_PER_TASK \
--pairtoolspath $PAIRTOOLS \
--juicerpath $JUICERPATH
```

2. Submit the script on biowulf: 
 
 
```
[biowulf]$ sbatch --mem=40g --cpus-per-task=8  myscript
```



Using Swarm

Using the 'swarm' utility, one can submit many jobs to the cluster to run concurrently. 


Set up a swarm command file (eg /data/$USER/cmdfile). 



```

cd /data/$USER/dir1; hint pre -d $TESTDIR/test.bam ...
cd /data/$USER/dir2; hint pre -d $TESTDIR/test.bam ...
cd /data/$USER/dir3; hint pre -d $TESTDIR/test.bam ...
...
cd /data/$USER/dir20; hint pre -d $TESTDIR/test.bam ...

```


 submit the swarm job:
 
```
$ swarm -f cmdfile --module hint -g 40 -t 16
```

For more information regarding running swarm, see <swarm.html>



















