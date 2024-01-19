

document.querySelector('title').textContent = 'xTea: comprehensive transposable element analyzer';
**xTea: comprehensive transposable element analyzer**


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



xTea (x-Transposable element analyzer), is a tool for identifying TE insertions
in whole-genome sequencing data. Whereas existing methods are mostly designed for shortread
data, xTea can be applied to both short-read and long-read data. xTea outperforms other short read-based methods 
for both germline and somatic TE insertion discovery. 



### References:


* Chong Chu, Rebeca Borges-Monroy, Vinayak V. Viswanadham, Soohyun Lee, Heng Li,
Eunjung Alice Lee and Peter J. Park,   

*Comprehensive identification of transposable element insertions using multiple sequencing technologies*   

[Nature Communications](https://www.nature.com/articles/s41467-021-24041-8) **12**, Article number: 3836 (2021).


Documentation
* [xTea Github page](https://github.com/parklab/xTea)
* [xTea Demo page](https://github.com/parklab/xTea/blob/master/Demo/demo_readme.md)


Important Notes
* Module Name: xTea (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set
	+ **XTEA\_HOME**  installation directory
	+ **XTEA\_BIN**       executable directory
	+ **XTEA\_DATA**  sample data directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=8g -c8 --gres=lscratch:10**
[user@cig 3335 ~]$ **module load xTea**
[+] Loading singularity  3.10.5  on cn3335
[+] Loading xTea  1.0.0 
[user@cn3335 ~]$ **xtea -h**
Usage: xtea [options]

Options:
  -h, --help            show this help message and exit
  -D, --decompress      Decompress the rep lib and reference file
  -M, --mosaic          Calling mosaic events from high coverage data
  -C, --case_control    Run in case control mode
  --denovo              Run in de novo mode
  -U, --user            Use user specific parameters instead of automatically
                        calculated ones
  --force               Force to start from the very beginning
  --hard                This is hard-cut for fitering out coverage abnormal
                        candidates
  --tumor               Working on tumor samples
  --purity=PURITY       Tumor purity
  --lsf                 Indiates submit to LSF system
  --slurm               Indiates submit to slurm system
  --resume              Resume the running, which will skip the step if output
                        file already exists!
  -V, --version         Print xTea version
  -i FILE, --id=FILE    sample id list file
  -a FILE, --par=FILE   parameter file
  -l FILE, --lib=FILE   TE lib config file
  -b FILE, --bam=FILE   Input bam file
  -x FILE, --x10=FILE   Input 10X bam file12878
  -n CORES, --cores=CORES
                        number of cores
  -m MEMORY, --memory=MEMORY
                        Memory limit in GB
  -q PARTITION, --partition=PARTITION
                        Which queue to run the job
  -t TIME, --time=TIME  Time limit
  -p WFOLDER, --path=WFOLDER
                        Working folder
  -r REF, --ref=REF     reference genome
  -g GENE, --gene=GENE  Gene annotation file
  --xtea=XTEA           xTEA folder
  -f FLAG, --flag=FLAG  Flag indicates which step to run (1-clip, 2-disc,
                        4-barcode, 8-xfilter, 16-filter, 32-asm)
  -y REP_TYPE, --reptype=REP_TYPE
                        Type of repeats working on: 1-L1, 2-Alu, 4-SVA,
                        8-HERV, 16-Mitochondrial
  --flklen=FLKLEN       flank region file
  --nclip=NCLIP         cutoff of minimum # of clipped reads
  --cr=CLIPREP          cutoff of minimum # of clipped reads whose mates map
                        in repetitive regions
  --nd=NDISC            cutoff of minimum # of discordant pair
  --nfclip=NFILTERCLIP  cutoff of minimum # of clipped reads in filtering step
  --nfdisc=NFILTERDISC  cutoff of minimum # of discordant pair of each sample
                        in filtering step
  --teilen=TEILEN       minimum length of the insertion for future analysis
  -o FILE, --output=FILE
                        The output file
  --blacklist=FILE      Reference panel database for filtering, or a blacklist
                        region
[user@cn3335 ~]$ **ls $XTEA\_BIN**
python  shell  xtea  xtea_hg19  xtea_long

```

Prepare the sample data data to be used:

```

[user@cn3335 ~]$ **ln -s $XTEA\_DATA/NA12878\_S1.bam**
[user@cn3335 ~]$ **samtools index NA12878\_S1.bam**

[user@cn3335 ~]$ **ln -s $XTEA\_DATA/gencode.v33.chr\_patch\_hapl\_scaff.basic.annotation.gff3**
[user@cn3335 ~]$ **ln -s /fdb/igenomes/Homo\_sapiens/Ensembl/GRCh37/Sequence/WholeGenomeFasta/genome.fa** 
[user@cn3335 ~]$ **ln -s /fdb/igenomes/Homo\_sapiens/Ensembl/GRCh37/Sequence/WholeGenomeFasta/genome.fa.index** 
[user@cn3335 ~]$ **ln -s /fdb/igenomes/Homo\_sapiens/Ensembl/GRCh37/Sequence/WholeGenomeFasta/genome.fa.fai** 
[user@cn3335 ~]$ **cp $XTEA\_DEMO/\* .**

```

Create a script to be submitted to the cluster:

```

[user@cn3335 ~]$ **sh run\_gnrt\_pipeline.sh**

```

Submit the script: 

```

[user@cn3335 ~]$ **source submit\_jobs.sh**
59538871
59538874
59538877

```





