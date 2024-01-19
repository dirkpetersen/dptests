

document.querySelector('title').textContent = 'kb-python: A wrapper for the kallisto | bustools workflow for single-cell RNA-seq pre-processing ';
**kb-python: A wrapper for the kallisto | bustools workflow for single-cell RNA-seq pre-processing**


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



kb-python is a python package for processing single-cell RNA-sequencing. It wraps the kallisto | bustools single-cell RNA-seq command line tools in order to unify multiple processing workflows.

- gget enables efficient querying of genomic reference databases to support the analysis of sequencing data. 
- ffq: A tool to find sequencing data and metadata from public databases. 
### References:


	
	Melsted, P., Booeshaghi, A.S., et al.   
	
	Modular, efficient and constant-memory single-cell RNA-seq preprocessing.  *Nat Biotechnol 39, 813–818 (2021).*    
	
	 [journal](https://doi.org/10.1038/s41587-021-00870-2)Documentation
	* [kb-python GitHub Page](https://github.com/pachterlab/kb_python)
	* [gget GitHub Page](https://github.com/pachterlab/gget)
	* [ffq GitHub Page](https://github.com/pachterlab/ffq)Important Notes
	* Module Name: kb-python (see [the modules page](/apps/modules.html) for more information)
	* gget/0.3.11 and ffq/0.3.0 are installed with the module.
Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=4g --gres=lscratch:10** 
[user@cn3144 ]$ **module load kb-python** 
[+] Loading singularity
[+] Loading kb-python  0.27.3
[user@cn3144 ]$ **kb**
usage: kb [-h] [--list]  ...

kb\_python 0.27.3

positional arguments:
 
 info Display package and citation information
 compile Compile `kallisto` and `bustools` binaries from source
 ref Build a kallisto index and transcript-to-gene mapping
 count Generate count matrices from a set of single-cell FASTQ files

options:
 -h, --help Show this help message and exit
 --list Display list of supported single-cell technologies

[user@cn3144 ]$ **kb info**
kb\_python 0.27.3
kallisto: 0.48.0 (/opt/conda/lib/python3.10/site-packages/kb\_python/bins/linux/kallisto/kallisto)
bustools: 0.41.0 (/opt/conda/lib/python3.10/site-packages/kb\_python/bins/linux/bustools/bustools)
kb is a python package for rapidly pre-processing single-cell RNA-seq data. It
is a wrapper for the methods described in [2].

The goal of the wrapper is to simplify downloading and running of the kallisto
[1] and bustools [2] programs. It was inspired by Sten Linnarsson’s loompy
fromfq command (http://linnarssonlab.org/loompy/kallisto/index.html)

The kb program consists of two parts:

The `kb ref` command builds or downloads a species-specific index for
pseudoalignment of reads. This command must be run prior to `kb count`, and it
runs the `kallisto index` [1].

The `kb count` command runs the kallisto [1] and bustools [2] programs. It can
be used for pre-processing of data from a variety of single-cell RNA-seq
technologies, and for a number of different workflows (e.g. production of gene
count matrices, RNA velocity analyses, etc.). The output can be saved in a
variety of formats including mix and loom. Examples are provided below.

Examples are available at: https://www.kallistobus.tools/tutorials

References
==========
[1] Bray, N. L., Pimentel, H., Melsted, P., & Pachter, L. (2016). Near-optimal
probabilistic RNA-seq quantification. Nature biotechnology, 34(5), 525.
[2] Melsted, P., Booeshaghi, A. S., Liu, L., Gao, F., Lu, L., Min, K. H., da
Veiga Beltrame, E., Hjorleifsson, K. E., Gehring, J., & Pachter, L. (2021).
Modular and efficient pre-processing of single-cell RNA-seq. Nature
Biotechnology.

[user@cn3144 ]$ **gget**
usage: gget [-h] [-v] {ref,search,info,seq,muscle,blast,blat,enrichr,archs4,setup,alphafold,pdb} ...

gget v0.3.11

positional arguments:
 {ref,search,info,seq,muscle,blast,blat,enrichr,archs4,setup,alphafold,pdb}
 ref Fetch FTPs for reference genomes and annotations by species.
 search Fetch gene and transcript IDs from Ensembl using free-form search terms.
 info Fetch gene and transcript metadata using Ensembl IDs.
 seq Fetch nucleotide or amino acid sequence (FASTA) of a gene (and all isoforms) or transcript by Ensembl, WormBase or FlyBase ID.
 muscle Align multiple nucleotide or amino acid sequences against each other (using the Muscle v5 algorithm).
 blast BLAST a nucleotide or amino acid sequence against any BLAST database.
 blat BLAT a nucleotide or amino acid sequence against any BLAT UCSC assembly.
 enrichr Perform an enrichment analysis on a list of genes using Enrichr.
 archs4 Find the most correlated genes or the tissue expression atlas of a gene using data from the human and mouse RNA-seq database
 ARCHS4 (https://maayanlab.cloud/archs4/).
 setup Install third-party dependencies for a specified gget module.
 alphafold Predicts the structure of a protein using a slightly simplified version of AlphaFold v2.1.0
 (https://doi.org/10.1038/s41586-021-03819-2).
 pdb Query RCSB PDB for the protein structutre/metadata of a given PDB ID.

options:
 -h, --help Print manual.
 -v, --version Print version.

[user@cn3144 ]$ **ffq**
usage: ffq [-h] [-o OUT] [-l LEVEL] [--ftp] [--aws] [--gcp] [--ncbi] [--split] [--verbose] [--version] IDs [IDs ...]

ffq 0.3.0: A command line tool to find sequencing data from SRA / GEO / ENCODE / ENA / EBI-EMBL / DDBJ / Biosample.

positional arguments:
 IDs One or multiple SRA / GEO / ENCODE / ENA / EBI-EMBL / DDBJ / Biosample accessions, DOIs, or paper titles

options:
 -h, --help Show this help message and exit
 -o OUT Path to write metadata (default: standard out)
 -l LEVEL Max depth to fetch data within accession tree
 --ftp Return FTP links
 --aws Return AWS links
 --gcp Return GCP links
 --ncbi Return NCBI links
 --split Split output into separate files by accession (`-o` is a directory)
 --verbose Print debugging information
 --version show program's version number and exit

[user@cn3144 ]$ **kb ref -i index.idx -g t2g.txt -f1 transcriptome.fa $(gget ref --ftp -w dna,gtf drosophila\_melanogaster)**
[user@cn3144 ]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```


Batch job
Most jobs should be run as [batch jobs](/docs/userguide.html#submit).
Create a batch input file (e.g. kb-python.sh) similar to the following example:



```

#!/bin/bash

module load kb-python
kb ref -i index.idx -g t2g.txt -f1 transcriptome.fa $(gget ref --ftp -w dna,gtf homo_sapiens)
kb count -i index.idx -g t2g.txt -x 10xv3 -o out $(ffq --ftp SRR10668798 | jq -r '.[] | .url' | tr '\n' ' ')

```

Submit this job using the Slurm [sbatch](/docs/userguide.html) command.



```
sbatch --cpus-per-task=6 --mem=4g kb-python.sh
```

Swarm of Jobs 
A [swarm of jobs](/apps/swarm.html) is an easy way to submit a set of independent commands requiring identical resource
s.
Create a swarmfile (e.g. kb-python.swarm). For example:



```

kb count -i index.idx -g t2g.txt -o out1/ -x 10xv3 read1_R1.fastq.gz read1_R2.fastq.gz
kb count -i index.idx -g t2g.txt -o out2/ -x 10xv3 read2_R1.fastq.gz read2_R2.fastq.gz

```

Submit this job using the [swarm](/apps/swarm.html) command.



```
swarm -f kb-python.swarm -g 4 -t 6 --module kb-python
```

where


|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| -g #  Number of Gigabytes of memory required for each process (1 line in the swarm command file)
 | -t #  Number of threads/CPUs required for each process (1 line in the swarm command file).
 | --module kb-python  Loads the kb-python module for each subjob in the swarm
 | |
 | |
 | |


