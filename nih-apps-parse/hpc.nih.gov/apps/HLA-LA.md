

document.querySelector('title').textContent = 'HLA-LA on Biowulf';
HLA-LA on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Additional References](#ref)
[Interactive job](#int) 
[Batch job](#sbatch) 
 |



HLA\*LA carries out HLA typing based on a population reference graph and employs a new linear projection method to align reads to the graph.

Previously called HLA\*PRG:LA, the application was developed by Alexander Dilthey at NHGRI.

### Reference:


* Alexander T Dilthey, Alexander J Mentzer, Raphael Carapito, Clare Cutland, Nezih Cereb, Shabir A Madhi, Arang Rhie, Sergey Koren, Seiamak Bahram, Gil McVean, Adam M Phillippy
 [**HLA\*LA—HLA typing from linearly projected graph alignments.**](https://academic.oup.com/bioinformatics/article/35/21/4394/5426702)*Bioinformatics, Volume 35, Issue 21, 1 November 2019, Pages 4394–4396*






Documentation
* [Github Site](https://github.com/DiltheyLab/HLA-LA/)
* [Description of the algorithm](https://genomeinformatics.github.io/HLA-PRG-LA/)


Important Notes
* Module Name: HLA-LA (see [the modules page](/apps/modules.html) for more information)
* Multithreaded.
* Environment variables set:
+ HLA\_LA\_GRAPHS HLA\*LA Graph References
+ HLA\_LA\_TESTDATA HLA\*LA Test Data


Additional References

Please contact the HPC staff (staff@hpc.nih.gov) if you want additional reference files installed for HLA-LA. 
Note that HLA-LA will break for all users if your reference file is incorrectly formatted.
See [this documentation](https://github.com/DiltheyLab/HLA-LA/#adding-further-references) for the format.
We highly recommend using an already installed reference as a template. For example: 



```

[user@cn3144 ~] module load HLA-LA
[user@cn3144 ~] cd $HLA_LA_GRAPHS/PRG_MHC_GRCh38_withIMGT/knownReferences
[user@cn3144 knownReferences] head PRG_MHC_GRCh38_withIMGT.txt > ~/testgraph.txt
[user@cn3144 knownReferences] tail PRG_MHC_GRCh38_withIMGT.txt >> ~/testgraph.txt

```


Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --gres=lscratch:50 --cpus-per-task=8 --mem=60g** 
salloc.exe: Pending job allocation 46116226
salloc.exe: job 46116226 queued and waiting for resources
salloc.exe: job 46116226 has been allocated resources
salloc.exe: Granted job allocation 46116226
salloc.exe: Waiting for resource configuration
salloc.exe: Nodes cn3144 are ready for job

[user@cn3144 ~]$ **module load HLA-LA**
[+] Loading singularity  3.10.5  
[+] Loading HLA-LA  1.0.3 

[user@cn3144 ~]$ **cd /lscratch/$SLURM\_JOBID**

[user@cn3144 46116226]$ **cp $HLA\_LA\_TESTDATA/NA12878.mini.cram .**

[user@cn3144 46116226]$ **samtools index NA12878.mini.cram**

[user@cn3144 46116226]$ **HLA-LA.pl --BAM NA12878.mini.cram \
 --graph PRG\_MHC\_GRCh38\_withIMGT --sampleID NA12878 \
 --maxThreads 7 --workingDir .**
HLA-LA.pl

Identified paths:
        samtools_bin: /opt/conda/envs/hla-la/bin/samtools
        bwa_bin: /opt/conda/envs/hla-la/bin/bwa
        java_bin: /opt/conda/envs/hla-la/bin/java
        picard_sam2fastq_bin: /opt/conda/envs/hla-la/bin/picard
        General working directory: /lscratch/4506949
        Sample-specific working directory: /lscratch/4506949/NA12878

Using /opt/conda/envs/hla-la/opt/hla-la/src/../graphs/PRG_MHC_GRCh38_withIMGT/knownReferences/1000G_B38.txt as reference file.
Extract reads from 534 regions...
Extract unmapped reads...
Merging...
Indexing...
Extract FASTQ...
        /opt/conda/envs/hla-la/bin/picard SamToFastq VALIDATION_STRINGENCY=LENIENT I=/lscratch/4506949/NA12878/extraction.bam F=/lscratch/4506949/NA12878/R_1.fastq F2=/lscratch/4506949/NA12878/R_2.fastq FU=/lscratch/4506949/NA12878/R_U.fastq 2>&1

Now executing:
../bin/HLA-LA --action HLA --maxThreads 7 --sampleID NA12878 --outputDirectory /lscratch/4506949/NA12878 --PRG_graph_dir /opt/conda/envs/hla-la/opt/hla-la/src/../graphs/PRG_MHC_GRCh38_withIMGT --FASTQU /lscratch/4506949/NA12878/R_U.fastq.splitLongReads --FASTQ1 /lscratch/4506949/NA12878/R_1.fastq --FASTQ2 /lscratch/4506949/NA12878/R_2.fastq --bwa_bin /opt/conda/envs/hla-la/bin/bwa --samtools_bin /opt/conda/envs/hla-la/bin/samtools --mapAgainstCompleteGenome 1 --longReads 0
Set maxThreads to 7

[...]

[user@cn3144 46116226]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. HLA.sh). For example:



```

#!/bin/bash
set -e
cd /lscratch/$SLURM_JOBID
module load HLA-LA

cp /data/$USER/myfile.cram .
samtools index myfile.cram

cpus=$(( SLURM_CPUS_PER_TASK - 1 ))
echo "Running on $cpus CPUs"
HLA-LA.pl --BAM myfile.cram --graph PRG_MHC_GRCh38_withIMGT --sampleID myfile --maxThreads $cpus --workingDir .

# copy output from /lscratch back to /data area
cp -r myfile/hla  /data/$USER/


```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=32 --mem=100g --gres=lscratch:100 --time=1-00:00:00 HLA.sh
```









