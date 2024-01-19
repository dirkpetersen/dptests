

document.querySelector('title').textContent = 'parsnp: efficient microbial core genome alignment and SNP detection ';
**parsnp: efficient microbial core genome alignment and SNP detection** 


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



Parsnp is a command-line-tool for efficient microbial core genome alignment and SNP detection. 
Parsnp was designed to work in tandem with Gingr, a flexible platform for visualizing genome alignments and phylogenetic trees



### References:


* Todd J Treangen, Brian D Ondov, Sergey Koren and Adam M Phillippy   

*The Harvest suite for rapid core-genome alignment and visualization of thousands of intraspecific microbial genomes*   

[Genome Biology](https://link.springer.com/article/10.1186/s13059-014-0524-x), 2014, 15:524. http://genomebiology.com/2014/15/11/524


Documentation
* [Parsnp Github page](https://github.com/marbl/parsnp)
* [Harvest Project page](https://harvest.readthedocs.io/en/latest/)


Important Notes
* Module Name: parsnp (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Multithreaded
* Unusual environment variables set
	+ **PARSNP\_HOME**  installation directory
	+ **PARSNP\_BIN**       executable directory
	+ **PARSNP\_DATA**  sample data directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=16g --cpus-per-task=16**
[user@cn3335 ~]$**module load parsnp** 
[+] Loading singularity  3.10.3  on cn3335
[+] Loading parsnp 1.7.4  ...
user@cn3335 ~]$ **parsnp -h** 
|--Parsnp v1.2--|
For detailed documentation please see --> http://harvest.readthedocs.org/en/latest
usage: parsnp [options] [-g|-r|-q](see below) -d <genome_dir> -p <threads>

Parsnp quick start for three example scenarios:
1) With reference & genbank file:
 >parsnp -g <reference_genbank_file1,reference_genbank_file2,..> -d <genome_dir> -p <threads>

2) With reference but without genbank file:
 >parsnp -r <reference_genome> -d <genome_dir> -p <threads>

3) Autorecruit reference to a draft assembly:
 >parsnp -q <draft_assembly> -d <genome_db> -p <threads>

[Input parameters]
<<input/output>>
 -c = <flag>: (c)urated genome directory, use all genomes in dir and ignore MUMi? (default = NO)
 -d = <path>: (d)irectory containing genomes/contigs/scaffolds
 -r = <path>: (r)eference genome (set to ! to pick random one from genome dir)
 -g = <string>: Gen(b)ank file(s) (gbk), comma separated list (default = None)
 -o = <string>: output directory? default [./P_CURRDATE_CURRTIME]
 -q = <path>: (optional) specify (assembled) query genome to use, in addition to genomes found in genome dir (default = NONE)

<<MUMi>>
 -U = <float>: max MUMi distance value for MUMi distribution
 -M = <flag>: calculate MUMi and exit? overrides all other choices! (default: NO)
 -i = <float>: max MUM(i) distance (default: autocutoff based on distribution of MUMi values)

<<MUM search>>
 -a = <int>: min (a)NCHOR length (default = 1.1*Log(S))
 -C = <int>: maximal cluster D value? (default=100)
 -z = <path>: min LCB si(z)e? (default = 25)

<<LCB alignment>>
 -D = <float>: maximal diagonal difference? Either percentage (e.g. 0.2) or bp (e.g. 100bp) (default = 0.12)
 -e = <flag> greedily extend LCBs? experimental! (default = NO)
 -n = <string>: alignment program (default: libMUSCLE)
 -u = <flag>: output unaligned regions? .unaligned (default: NO)

<<Recombination filtration>>
 -x = <flag>: enable filtering of SNPs located in PhiPack identified regions of recombination? (default: NO)

<<Misc>>
 -h = <flag>: (h)elp: print this message and exit
 -p = <int>: number of threads to use? (default= 1)
 -P = <int>: max partition size? limits memory usage (default= 15000000)
 -v = <flag>: (v)erbose output? (default = NO)
 -V = <flag>: output (V)ersion and exit

```

Download sample data to the current folder:

```

[user@cn3335 ~]$**cp -r $PARSNP\_DATA/\* .**

```

Run parsnp on the sample data: 

```

[user@cn3335 ~]$ **parsnp -g ref/England1.gbk -d genomes -c**
|--Parsnp v1.2--|
For detailed documentation please see --> http://harvest.readthedocs.org/en/latest
*****************************************************************************
SETTINGS:
|-refgenome:    ref/England1.gbk.fna
|-aligner:      libMUSCLE
|-seqdir:       genomes
|-outdir:       /data/user/parsnp/P_2022_10_20_091716517037
|-OS:           Linux
|-threads:      32
*****************************************************************************

<<Parsnp started>>

-->Reading Genome (asm, fasta) files from genomes..
  |->[OK]
-->Reading Genbank file(s) for reference (.gbk) ref/England1.gbk..
  |->[OK]
-->Running Parsnp multi-MUM search and libMUSCLE aligner..
  |->[OK]
-->Running PhiPack on LCBs to detect recombination..
  |->[SKIP]
-->Reconstructing core genome phylogeny..
  |->[OK]
-->Creating Gingr input file..
  |->[OK]
-->Calculating wall clock time..
  |->Aligned 47 genomes in 0.65 seconds

<<Parsnp finished! All output available in /data/user/parsnp/P_2022_10_20_091716517037>>

Validating output directory contents...
        1)parsnp.tree:          newick format tree                      [OK]
        2)parsnp.ggr:           harvest input file for gingr (GUI)      [OK]
        3)parsnp.xmfa:          XMFA formatted multi-alignment          [OK]
[user@cn3335 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```





