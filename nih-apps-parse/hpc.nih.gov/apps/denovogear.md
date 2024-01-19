

document.querySelector('title').textContent = 'DeNovoGear: Estimating de novo mutations from related individuals and cells ';
**DeNovoGear: Estimating de novo mutations from related individuals and cells** 


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



DeNovoGear is az software for analyzing de novo mutations from familial and somatic tissue sequencing data. 
It uses likelihood-based error modeling to reduce the false positive rate of mutation discovery in exome analysis 
and fragment information to identify the parental origin of germ-line mutations. DeNovoGear has been used on human
whole-genome sequencing data to produce a set of predicted de novo insertion and/or deletion (indel) mutations.



### Reference:


* Avinash Ramu, Michiel J Noordam, Rachel S Schwartz, Arthur Wuster, Matthew E Hurles, Reed A Cartwright and Donald F Conrad  

*DeNovoGear: de novo indel and point mutation discovery and phasing*   

[Nature Methods](https://www.nature.com/articles/nmeth.2611), 2013, **10**(10), p.985.

Documentation
	+ [DeNovoGear Github page](https://github.com/ultimatesource/denovogear)Important Notes
	+ Module Name: denovogear (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
	+ Unusual environment variables set
		- **DENOVOGEAR\_HOME**  installation directory
		- **DENOVOGEAR\_BIN**       executable directory
		- **DENOVOGEAR\_SRC**       source code directory
		- **DENOVOGEAR\_DATA**  sample data directory
Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=12g -c5 --gres=lscratch:10**
[user@cn3335 ~]$ **module load denovogear** 
[+] Loading samtools 1.17  ...
[+] Loading denovogear 1.1.1  ...
[user@cn3335 ~]$
[user@cn3335 ~]$ **dng**
USAGE: dng command [options]
       dng help
       dng help [command]
[user@cn3335 ~]$ **dng-dnm**
DeNovoGear v1.1.1

Usage:
Autosomes:
        dng dnm auto --bcf bcf_f --ped ped_f [OR] dng dnm auto --vcf vcf_f --ped ped_f
X chromosome in male offspring:
        dng dnm XS --bcf bcf_f --ped ped_f [OR] dng dnm XS --vcf vcf_f --ped ped_f
X chromosome in female offspring:
        dng dnm XD --bcf bcf_f --ped ped_f [OR] dng dnm XD --vcf vcf_f --ped ped_f

Input:
DNM:
--ped:   Ped file to describe relationship between the samples.
--bcf:   BCF file, contains per-sample read depths and genotype likelihoods.
--vcf:   VCF file, contains per-sample read depths and genotype likelihoods.
Phaser:
--dnm: Tab delimited list of denovo mutations to be phased, format: chr pos inherited_base denovo_base.[example: 1 2000 A C]
--pgt: Tab delimited genotypes of child and parents at SNP sites near denovo sites, format: chr pos GT_child GT_parent1 GT_parent2.[example: 1 2000 AC AC AA]
--bam: alignment file (.bam) of the child.
--window: optional argument which is the maximum distance between the DNM and a phasing site. The default value is 1000.

Output:
--output_vcf:    vcf file to write the output to.

Parameters:
--snp_mrate:     Mutation rate prior for SNPs. [1e-8]
--indel_mrate:   Mutation rate prior for INDELs. [1e-9]
--pair_mrate:    Mutation rate prior for paired sample analysis. [1e-9]
--indel_mu_scale:        Scaling factor for indel mutation rate. [1]
--pp_cutoff:     Posterior probability threshold. [0.0001]
--rd_cutoff:     Read depth filter, sites where either one of the sample have read depth less than this threshold are filtered out. [10]
--region:        Region of the BCF file to perform denovo calling. [string of the form "chr:start-end"
[user@cn3335 ~]$ **dng-call**
Usage:
  dng call [options] input1 input2 input3 ...

Allowed Options:
  -f [ --fasta ] arg                 faidx indexed reference sequence file
  -l [ --min-qlen ] arg (=0)         minimum query length
  -m [ --min-prob ] arg (=0.1)       minimum probability for reporting a
                                     mutation
  --mu arg (=1e-9)                   the germline mutation rate
  --mu-somatic arg (=0)              the somatic mutation rate
  --mu-library arg (=0)              the library prep mutation rate
  --nuc-freqs arg (=0.3,0.2,0.2,0.3) nucleotide frequencies in ACGT order
  -p [ --ped ] arg                   the pedigree file
  -q [ --min-basequal ] arg (=0)     minimum base quality
  -Q [ --min-mapqual ] arg (=0)      minimum mapping quality
  -r [ --region ] arg                chromosomal region
  -R [ --ref-weight ] arg (=1)       weight given to reference base for
                                     population prior
  -s [ --sam-files ] arg             file containing a list of input filenames,
                                     one per line
  --theta arg (=0.001)               the population diversity
  -o [ --output ] arg (=-)           Output VCF/BCF file
  --version                          display version information
  --help                             display usage informaiton
  --arg-file arg                     read command-line arguments from a file
[user@cn3335 ~]$ 


```

etc.   
   

End the interactive session:

```

[user@cn3335 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```
