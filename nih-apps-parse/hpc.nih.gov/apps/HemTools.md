

document.querySelector('title').textContent = 'HemTools: a collection of NGS pipelines and bioinformatic analyses.';
**HemTools: a collection of NGS pipelines and bioinformatic analyses.**


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



HemTools is a collection of NGS pipelines and bioinformatic analysis tools. It includes tools for data visualization, motif analysis, integrative analysis, bioinformatica analysis, differential analysis, CRISPR analysis, and more.



### References:


* Phillip A. Doerfler, Ruopeng Feng, Yichao Li, Lance E. Palmer, Shaina N. Porter, Henry W. Bell, Merlin Crossley, Shondra M. Pruett-Miller, Yong Cheng and Mitchell J. Weiss,   

*Activation of γ-globin gene expression by GATA1 and NF-Y in hereditary persistence of fetal hemoglobin*   

[Nature Genetics](https://www.nature.com/articles/s41588-021-00904-0) 2021, vol.53, p.1177–1186.


Documentation
* [HemTools Github page](https://github.com/YichaoOU/HemTools)
* [HemTools Documentation](https://hemtools.readthedocs.io/en/latest/)


Important Notes
* Module Name: HemTools (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set
	+ **HEMTOOLS\_HOME**  installation directory
	+ **HEMTOOLS\_BIN**  executable directory
	+ **HEMTOOLS\_DIR**  source code directory
	+ **HEMTOOLS\_DATA**  configuration files and models directory
	+ **HEMTOOLS\_BENCHMARKS**  benchmark datasets directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=16g --gres=gpu:p100:1,lscratch:10 -c4** 
[user@cn3104 ~]$ **module load HemTools**
[+] Loading singularity  3.10.5  on cn3104
[+] Loading HemTools 20210512
[user@cn3104 ~]$ **HemTools -h**
usage: HemTools [-h] [-v]
                {cut_run,cut_run_histone,chip_seq_pair,chip_seq_single,atac_seq,report_bug,rna_seq,my_dir,volcano_plot,crispr_seq}
                ...

HemTools: performs NGS pipelines and other common analyses. Contact:
Yichao.Li@stjude.org or Yong.Cheng@stjude.org

positional arguments:
  {cut_run,cut_run_histone,chip_seq_pair,chip_seq_single,atac_seq,report_bug,rna_seq,my_dir,volcano_plot,crispr_seq}
                        Available APIs in HemTools
    cut_run             CUT & RUN pipeline
    cut_run_histone     CUT & RUN pipeline
    chip_seq_pair       Paired-end ChIP-seq pipeline
    chip_seq_single     Single-end ChIP-seq pipeline
    atac_seq            ATAC-seq pipeline
    report_bug          Email the log files to the developer.
    rna_seq             RNA-seq pipeline
    my_dir              CD, search, and list my dirs
    volcano_plot        Data visualization: Volcano plot
    crispr_seq          Genome-wide CRISPR Screening pipeline
[user@cn3104 ~]$ **HemTools cut\_run -h**
usage: HemTools cut_run [-h] [-j JID] [--short] [--debug] [-f INPUT]
                        [-d DESIGN_MATRIX] [--guess_input] [-i INDEX_FILE]
                        [-g GENOME] [-b BLACKLIST] [-s CHROM_SIZE]
                        [-e EFFECTIVEGENOMESIZE]

optional arguments:
  -h, --help            show this help message and exit
  -j JID, --jid JID     enter a job ID, which is used to make a new directory.
                        Every output will be moved into this folder. (default:
                        {{subcmd}}_$USER_2023-05-30)
  --short               Force to use the short queue. (only if R1+R2 fastq.gz
                        size <=250M) (default: False)
  --debug               Not for end-user. (default: False)
  -f INPUT, --input INPUT
                        tab delimited 3 columns (tsv file): Read 1 fastq, Read
                        2 fastq, sample ID (default: None)
  -d DESIGN_MATRIX, --design_matrix DESIGN_MATRIX
                        tab delimited 3 columns (tsv file): treatment sample
                        ID, control sample ID, peakcall ID (default: None)
  --guess_input         Let the program generate the input files for you.
                        (default: False)

Genome Info:
  -i INDEX_FILE, --index_file INDEX_FILE
                        BWA index file (default:
                        /opt/conda/envs/hemtools/HemTools-
                        HemTools/subcmd/../hg19/bwa_16a_index/hg19.fa)
  -g GENOME, --genome GENOME
                        genome version: hg19, hg38, mm10, mm9. (default: hg19)
  -b BLACKLIST, --Blacklist BLACKLIST
                        Blacklist file (default:
                        /opt/conda/envs/hemtools/HemTools-
                        HemTools/subcmd/../hg19/Hg19_Blacklist.bed)
  -s CHROM_SIZE, --chrom_size CHROM_SIZE
                        chrome size (default:
                        /opt/conda/envs/hemtools/HemTools-
                        HemTools/subcmd/../hg19/hg19.chrom.sizes)
  -e EFFECTIVEGENOMESIZE, --effectiveGenomeSize EFFECTIVEGENOMESIZE
                        effectiveGenomeSize for bamCoverage (default:
                        2451960000)

```

Exit the application:   


```

[user@cn3104 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```





